#!/usr/bin/env python3
"""Generate the committed Fortran sources from the fypp templates.

    python3 scripts/fypp_deploy.py                 # regenerate src/ and test/
    python3 scripts/fypp_deploy.py --check         # fail if the committed tree is out of date

Every fypp/src/*.fypp becomes src/<stem>.f90 (or .F90 when the template carries cpp
directives) and every fypp/test/**/*.fypp becomes test/<same relative path>.f90.  fypp runs
with -I include, the result goes through fprettify with the project flags and through the
labelled-continue realignment that undoes fprettify's re-indentation of numbered statements.
"""

import argparse
import concurrent.futures
import difflib
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

FPRETTIFY_FLAGS = [
    "-i", "4", "-l", "132", "-w", "2", "--disable-indent", "--strip-comments",
    "--c-relations", "--enable-replacements", "--enable-decl", "--whitespace-comma", "0",
]

# Templates that the pipeline does not own yet.
NOT_YET_TEMPLATED = {}

# Templates whose committed output was edited by hand afterwards: regenerating them would revert
# those edits, so they need a reconciliation change of their own before they rejoin the pipeline.
SOURCE_DIVERGED = {}

EXCLUDED = dict(NOT_YET_TEMPLATED, **SOURCE_DIVERGED)

# Committed names that do not follow the extension rule below.
NAME_OVERRIDE = {}

CPP_DIRECTIVE = re.compile(r"^#\s*(if|ifdef|ifndef|elif|else|endif|define|undef)\b", re.M)
LABELLED_CONTINUE = re.compile(r"^\s*(\d+)\s+(continue)\s*$")


def output_name(stem, text):
    """Committed file name for a template: .F90 when it holds cpp directives, .f90 otherwise."""
    if stem in NAME_OVERRIDE:
        return NAME_OVERRIDE[stem]
    return stem + (".F90" if CPP_DIRECTIVE.search(text) else ".f90")


def realign_labelled_continue(path):
    """Restore the indentation fprettify strips from numbered continue statements."""
    with open(path) as fid:
        lines = [line.rstrip("\n") for line in fid]
    changed = False
    for i, line in enumerate(lines):
        match = LABELLED_CONTINUE.match(line)
        if not match or i == 0:
            continue
        previous = lines[i - 1]
        indent = len(previous) - len(previous.lstrip())
        new = " " * indent + match.group(1) + " " + match.group(2)
        if new != line:
            lines[i] = new
            changed = True
    if changed:
        with open(path, "w") as fid:
            fid.write("\n".join(lines) + "\n")


def _case_clash(dest_rel):
    """The committed path that differs from `dest_rel` in case only, on any filesystem."""
    folder = os.path.join(ROOT, os.path.dirname(dest_rel))
    want = os.path.basename(dest_rel)
    for name in os.listdir(folder):
        if name != want and name.lower() == want.lower():
            return os.path.join(os.path.dirname(dest_rel), name)
    return None


def generate(source, dest_dir, defines):
    """Run fypp + fprettify for one template; return (relative destination path, error)."""
    stem = os.path.splitext(os.path.basename(source))[0]
    with open(source) as fid:
        text = fid.read()
    if source.startswith(os.path.join(ROOT, "fypp", "test")):
        rel = os.path.relpath(source, os.path.join(ROOT, "fypp", "test"))
        dest_rel = os.path.join("test", os.path.splitext(rel)[0] + ".f90")
    else:
        dest_rel = os.path.join("src", output_name(stem, text))
        clash = _case_clash(dest_rel)
        if clash:
            return dest_rel, ("committed as %s; add the name to NAME_OVERRIDE" % clash)
    dest = os.path.join(dest_dir, dest_rel)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    cmd = ["fypp", "-I", os.path.join(ROOT, "include")]
    for name, value in defines:
        cmd += ["-D", "%s=%s" % (name, value)] if value is not None else ["-D", name]
    cmd += [source, dest]
    run = subprocess.run(cmd, capture_output=True, text=True)
    if run.returncode != 0:
        return dest_rel, "fypp: " + run.stderr.strip()
    run = subprocess.run(["fprettify"] + FPRETTIFY_FLAGS + [dest], capture_output=True, text=True)
    if run.returncode != 0:
        return dest_rel, "fprettify: " + run.stderr.strip()
    realign_labelled_continue(dest)
    return dest_rel, None


def collect_sources(only):
    sources, skipped = [], []
    src_dir = os.path.join(ROOT, "fypp", "src")
    for name in sorted(os.listdir(src_dir)):
        if not name.endswith(".fypp"):
            continue
        stem = name[:-len(".fypp")]
        if only and not any(re.fullmatch(pattern, stem) for pattern in only):
            continue
        if stem in EXCLUDED and not only:
            skipped.append((stem, EXCLUDED[stem]))
            continue
        sources.append(os.path.join(src_dir, name))
    test_dir = os.path.join(ROOT, "fypp", "test")
    if not only:
        for dirpath, _, names in sorted(os.walk(test_dir)):
            for name in sorted(names):
                if name.endswith(".fypp"):
                    sources.append(os.path.join(dirpath, name))
    return sources, skipped


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true",
                        help="do not write; fail if the committed tree differs")
    parser.add_argument("--jobs", type=int, default=os.cpu_count(),
                        help="templates processed in parallel")
    parser.add_argument("-D", dest="defines", action="append", default=[], metavar="NAME=VALUE",
                        help="passed through to fypp")
    parser.add_argument("--report", default=os.path.join(ROOT, "build", "fypp_deploy_check.diff"),
                        help="where --check writes the unified diff")
    parser.add_argument("--only", action="append", default=[], metavar="REGEX",
                        help="process only templates whose stem matches; bypasses the exclusions")
    args = parser.parse_args()

    defines = []
    for item in args.defines:
        name, _, value = item.partition("=")
        defines.append((name, value if _ else None))

    sources, skipped = collect_sources(args.only)
    for stem, reason in skipped:
        print("skip  %-20s %s" % (stem, reason))
    if not sources:
        print("no templates selected", file=sys.stderr)
        return 1

    staging = tempfile.mkdtemp(prefix="fypp_deploy_") if args.check else ROOT
    try:
        failures = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(generate, s, staging, defines): s for s in sources}
            produced = []
            for future in concurrent.futures.as_completed(futures):
                dest_rel, error = future.result()
                if error:
                    failures.append((dest_rel, error))
                else:
                    produced.append(dest_rel)
        if failures:
            for dest_rel, error in sorted(failures):
                print("FAIL  %s: %s" % (dest_rel, error), file=sys.stderr)
            return 1

        if not args.check:
            print("generated %d files" % len(produced))
            return 0

        differing, diff_text = [], []
        for dest_rel in sorted(produced):
            fresh = os.path.join(staging, dest_rel)
            committed = os.path.join(ROOT, dest_rel)
            new = open(fresh).read()
            old = open(committed).read() if os.path.exists(committed) else ""
            if new == old:
                continue
            differing.append(dest_rel)
            diff_text += list(difflib.unified_diff(
                old.splitlines(True), new.splitlines(True),
                fromfile="committed/" + dest_rel, tofile="generated/" + dest_rel))
        print("checked %d files, %d differ" % (len(produced), len(differing)))
        if differing:
            os.makedirs(os.path.dirname(args.report), exist_ok=True)
            with open(args.report, "w") as fid:
                fid.writelines(diff_text)
            for dest_rel in differing:
                print("OUT OF DATE  " + dest_rel, file=sys.stderr)
            print("unified diff: " + args.report, file=sys.stderr)
            return 1
        return 0
    finally:
        if args.check:
            shutil.rmtree(staging, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())

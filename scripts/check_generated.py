#!/usr/bin/env python3
"""Per-routine identity gate for the BLAS/LAPACK templating campaign.

    python3 scripts/check_generated.py --baseline origin/main --allow scripts/la_renames.tsv

Reads every la_* routine out of the baseline tree at a git ref and out of the working tree, then
asserts, name by name, that the working tree still holds the same code:

    PASS-BYTE   identical after collapsing whitespace
    PASS-NORM   identical after kind-role normalization and after undoing identifier renames
                that differ only by a leading kind letter (this is what accepts the tag and
                local-name repairs the templates make in the q and w instances)
    FAIL        anything else

Either side may carry exactly one `use la_constants_<kind>` line per routine, whose kind must
match the routine's own precision; the gate removes and verifies those lines rather than ignoring
them, on the baseline side too, so that a baseline which is itself already templated compares.
Names that legitimately change are listed in the allow-list file, one
`old_name  new_name|REMOVED|REFORMATTED|DOCTEXT  reason` row each; the last two name a body
difference of one routine, named as the working tree spells it, and the comparison the gate
holds it to.  The umbrella src/la_blas.F90 is compared as a
whole file, with only the `use` block and allow-listed procedure names permitted to differ.
"""

import argparse
import collections
import difflib
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import la_kindmap as K

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# The high-level API sits above the generated tiers and is not part of this comparison.
HIGH_LEVEL = {
    "la_cholesky", "la_constants", "la_determinant", "la_eigs", "la_eye", "la_inverse",
    "la_least_squares", "la_norms", "la_pinv", "la_qr", "la_schur", "la_solve", "la_state",
    "la_svd", "linear_algebra",
}
UMBRELLAS = {"src/la_blas.F90", "src/la_lapack.f90"}
USE_CONSTANTS = re.compile(r"(?m)^[ \t]*use la_constants_(sp|dp|qp)\b[^\n]*\n")
KIND_OF = {"s": "sp", "d": "dp", "q": "qp", "c": "sp", "z": "dp", "w": "qp"}


def git_show(ref, path):
    run = subprocess.run(["git", "-C", ROOT, "show", "%s:%s" % (ref, path)],
                         capture_output=True, text=True)
    return run.stdout if run.returncode == 0 else None


def git_files(ref):
    run = subprocess.run(["git", "-C", ROOT, "ls-tree", "-r", "--name-only", ref, "src/"],
                         capture_output=True, text=True, check=True)
    return [p for p in run.stdout.split("\n") if p]


def is_generated(path):
    name = os.path.basename(path)
    stem = os.path.splitext(name)[0]
    return (name.endswith((".f90", ".F90")) and stem.startswith("la_")
            and stem not in HIGH_LEVEL)


def baseline_tree(ref):
    out = {}
    for path in git_files(ref):
        if is_generated(path):
            out[path] = git_show(ref, path)
    return out


def working_tree():
    out = {}
    src = os.path.join(ROOT, "src")
    for name in sorted(os.listdir(src)):
        path = "src/" + name
        if is_generated(path):
            with open(os.path.join(src, name)) as fid:
                out[path] = fid.read()
    return out


def collect(tree):
    """`routine name -> (body, file)` over a whole tree, and the duplicate names found."""
    routines, duplicates = {}, []
    for path, text in sorted(tree.items()):
        if path in UMBRELLAS:
            continue
        for name, (chunk, _) in K.split_routines(text).items():
            if name in routines:
                duplicates.append((name, routines[name][1], path))
            routines[name] = (chunk, path)
    return routines, duplicates


def read_allow(path):
    allow = {}
    if not path:
        return allow
    with open(path) as fid:
        header = True
        for line in fid:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if header:
                header = False
                if line.split("\t")[0] == "old_name":
                    continue
            old, new, reason = (line.split("\t") + ["", ""])[:3]
            allow[old] = (new, reason)
    return allow


def strip_use(body, name):
    """Remove the imports of the per-kind constants and check they match the routine's kind."""
    found = USE_CONSTANTS.findall(body)
    if not found:
        return body, None
    expected = KIND_OF.get(_own_letter(name))
    for kind in found:
        if kind != expected:
            return USE_CONSTANTS.sub("", body), \
                "imports la_constants_%s, expected la_constants_%s" % (kind, expected)
    if len(found) > 1:
        return USE_CONSTANTS.sub("", body), "imports la_constants_%s %d times" % (found[0],
                                                                                  len(found))
    return USE_CONSTANTS.sub("", body), None


LETTER_PATTERNS = (r"(?:selctg|select)_([sdqczw])", r"i([sdqczw])(?:amax|max1)",
                   r"ila([sdqczw])(?:lc|lr|iag)", r"([sdqczw]).*")


def _own_letter(name):
    for pattern in LETTER_PATTERNS:
        m = re.fullmatch(pattern, name)
        if m:
            return m.group(1)
    return None


UPPER_TAG = re.compile(r"(?:@[RC][2LU]*I@)+(?=[A-Z])")
MARK = re.compile(r"@[A-Z0-9]+@")
QUOTED = re.compile(r"'[^'\n]*'")
LETTER_PAIRS = {(x, y): "@N@" for x in K.LETTERS for y in K.LETTERS if x != y}


def comment_at(line):
    """Index of the `!` that starts a comment, ignoring one inside a character literal."""
    quoted = False
    for i, ch in enumerate(line):
        if ch == "'":
            quoted = not quoted
        elif ch == "!" and not quoted:
            return i
    return len(line)


def mask_text(text):
    """Drop the kind role of every placeholder that sits in a comment or a character literal."""
    out = []
    for line in text.split("\n"):
        cut = comment_at(line)
        code = QUOTED.sub(lambda m: MARK.sub("@TAG@", m.group(0)), line[:cut])
        out.append(code + MARK.sub("@TAG@", line[cut:]))
    return "\n".join(out)


def no_space(text):
    """Every blank outside a character literal removed."""
    out, quoted = [], False
    for ch in text:
        if ch == "'":
            quoted = not quoted
        if quoted or not ch.isspace():
            out.append(ch)
    return "".join(out)


def no_comment(text):
    return "\n".join(line[:comment_at(line)] for line in text.split("\n"))


# Body differences accepted for one named routine each, with the comparison they are held to.
BODY_ALLOW = {
    "REFORMATTED": lambda a, b: no_space(a) == no_space(b),
    "DOCTEXT": lambda a, b: K.ws_norm(no_comment(a)) == K.ws_norm(no_comment(b)),
}


def compare(old_body, new_body, old_name, new_name, upper_bases, allowance=None):
    """Classify one routine and say which rung of the ladder accepted it.

    BYTE       the two bodies collapse to the same text
    NORM/kind  they agree once every precision token is replaced by its role
    NORM/name  they also need identifiers that differ only by a leading kind letter to be
               equated: the local names the q and w copies inherited from d and z
    NORM/tag   they also need the kind role of a placeholder that sits in a comment or in a
               character literal to be ignored: the stale 'DGETRF'-style tags and the stale kind
               words the q and w copies inherited.  Code outside a literal is never touched.
    ALLOW/...  the routine has a row in the allow-list naming a body difference and the
               comparison that row asks for accepts it.
    """
    if K.ws_norm(old_body) == K.ws_norm(new_body):
        return "BYTE", {}, None
    old_letter, new_letter = _own_letter(old_name), _own_letter(new_name)
    if old_letter is None or new_letter is None:
        if allowance and BODY_ALLOW[allowance](old_body, new_body):
            return "ALLOW/" + allowance.lower(), {}, None
        return "FAIL", {}, diff(old_body, new_body, old_name, new_name)
    a = K.normalize(old_body, old_letter, upper_bases)
    b = K.normalize(new_body, new_letter, upper_bases)
    if K.ws_norm(a) == K.ws_norm(b):
        return "NORM/kind", {}, None
    renames = K.letter_renames(a, b, LETTER_PAIRS)
    a2 = K.apply_renames(a, renames)
    b2 = K.apply_renames(b, K.letter_renames(b, a, LETTER_PAIRS))
    if K.ws_norm(a2) == K.ws_norm(b2):
        return "NORM/name", renames, None
    a3, b3 = mask_text(UPPER_TAG.sub("@TAG@", a2)), mask_text(UPPER_TAG.sub("@TAG@", b2))
    if K.ws_norm(a3) == K.ws_norm(b3):
        return "NORM/tag", renames, None
    if allowance and BODY_ALLOW[allowance](a3, b3):
        return "ALLOW/" + allowance.lower(), renames, None
    return "FAIL", {}, diff(K.ws_norm(a3), K.ws_norm(b3), old_name, new_name)


def diff(a, b, na, nb):
    return "\n".join(difflib.unified_diff(a.split("\n"), b.split("\n"),
                                          fromfile="baseline/" + na, tofile="working/" + nb,
                                          lineterm=""))


def compare_umbrella(old, new, allow, report):
    """The umbrella may differ only in its `use` block and in allow-listed procedure names."""
    use_block = re.compile(r"(?ms)^(     use .*?\n)(?=     implicit none)")
    old_use = use_block.search(old)
    new_use = use_block.search(new)
    if not old_use or not new_use:
        report.append("umbrella: could not locate the use block")
        return False
    old_rest = old[:old_use.start()] + old[old_use.end():]
    new_rest = new[:new_use.start()] + new[new_use.end():]
    for old_name, (new_name, _) in allow.items():
        if new_name not in BODY_ALLOW and new_name != "REMOVED":
            old_rest = re.sub(r"\b%s\b" % re.escape(old_name), new_name, old_rest)
    if old_use.group(1) != new_use.group(1):
        report.append("umbrella use block:\n" + diff(old_use.group(1), new_use.group(1),
                                                     "la_blas.F90", "la_blas.F90"))
    if old_rest == new_rest:
        return True
    report.append("umbrella body:\n" + diff(old_rest, new_rest, "la_blas.F90", "la_blas.F90"))
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", required=True, help="git ref holding the pre-refactor tree")
    parser.add_argument("--allow", default=None, help="tsv of intentional name changes")
    parser.add_argument("--report", default=os.path.join(ROOT, "build", "check_generated.txt"))
    args = parser.parse_args()

    allow = read_allow(args.allow)
    old_tree, new_tree = baseline_tree(args.baseline), working_tree()
    old_routines, old_dup = collect(old_tree)
    new_routines, new_dup = collect(new_tree)
    upper_bases = set(old_routines) | set(new_routines)

    report, counts = [], collections.Counter()
    failures = 0

    for name, first, second in old_dup:
        report.append("duplicate in the baseline: %s in %s and %s" % (name, first, second))
    for name, first, second in new_dup:
        report.append("duplicate in the working tree: %s in %s and %s" % (name, first, second))
        failures += 1

    bodies = {old[3:]: new for old, (new, _) in allow.items() if new in BODY_ALLOW}
    renames = {old: new for old, (new, _) in allow.items()
               if new != "REMOVED" and new not in BODY_ALLOW}
    removed = {old for old, (new, _) in allow.items() if new == "REMOVED"}
    matched_new = set()

    for name in sorted(old_routines):
        target = renames.get("la_" + name, "la_" + name)[3:]
        if "la_" + name in removed:
            counts["allow-removed"] += 1
            if name in new_routines:
                report.append("%s is allow-listed as REMOVED but still present" % name)
                failures += 1
            continue
        if target not in new_routines:
            report.append("missing: %s (expected as %s)" % (name, target))
            counts["missing"] += 1
            failures += 1
            continue
        matched_new.add(target)
        old_body, problem = strip_use(old_routines[name][0], name)
        if problem:
            report.append("baseline %s: %s" % (name, problem))
            failures += 1
        for was, now in renames.items():
            if re.search(r"\b%s\b" % re.escape(was), old_body):
                old_body = re.sub(r"\b%s\b" % re.escape(was), now, old_body)
                counts["call sites renamed"] += 1
        new_body, problem = strip_use(new_routines[target][0], target)
        if problem:
            report.append("%s: %s" % (target, problem))
            failures += 1
        verdict, applied, text = compare(old_body, new_body, name, target, upper_bases,
                                         bodies.get(target))
        key = "FAIL" if verdict == "FAIL" else \
            ("PASS-BYTE" if verdict == "BYTE" else "PASS-" + verdict)
        if verdict.startswith("ALLOW/"):
            report.append("%s %s: %s" % (key, target, allow["la_" + target][1]))
        counts[key] += 1
        if target != name:
            counts["allow-listed rename"] += 1
        if applied:
            report.append("%s %s: identifiers %s" % (key, target, ",".join(sorted(applied))))
        if verdict == "FAIL":
            failures += 1
            report.append("FAIL %s\n%s" % (target, text))

    for name in sorted(set(new_routines) - matched_new):
        report.append("extra: %s in %s is absent from the baseline" % (name, new_routines[name][1]))
        counts["extra"] += 1
        failures += 1

    for path in sorted(UMBRELLAS & set(old_tree)):
        if path not in new_tree:
            report.append("umbrella %s disappeared" % path)
            failures += 1
            continue
        if old_tree[path] == new_tree[path]:
            counts["umbrella PASS-BYTE"] += 1
        elif compare_umbrella(old_tree[path], new_tree[path], allow, report):
            counts["umbrella PASS-BYTE (use block and renames aside)"] += 1
        else:
            counts["umbrella FAIL"] += 1
            failures += 1

    width = max(len(k) for k in counts) if counts else 10
    print("baseline %s: %d routines in %d files" % (args.baseline, len(old_routines),
                                                    len(old_tree)))
    print("working tree: %d routines in %d files" % (len(new_routines), len(new_tree)))
    for key in sorted(counts):
        print("  %-*s %5d" % (width, key, counts[key]))
    if report:
        os.makedirs(os.path.dirname(args.report), exist_ok=True)
        with open(args.report, "w") as fid:
            fid.write("\n".join(report) + "\n")
        print("details: " + args.report)
    print("FAILURES: %d" % failures)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

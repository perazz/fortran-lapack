#!/usr/bin/env python3
"""Fence the optional precisions of every template behind their cpp guard.

    python3 scripts/guard_kinds.py            # rewrite the templates and the interface tables
    python3 scripts/guard_kinds.py --check    # fail if a template is not fenced

Two transformations, both idempotent, so the script can be re-run on a tree that already
carries them and on a tree that has moved on:

  templates   every outermost kind loop of fypp/src and fypp/test opens with the cpp fence
              LA_GUARD names for the kind of that iteration and closes with its #endif, so a
              qp instance is wrapped in #ifdef LA_WITH_QP and an xdp instance in
              #ifdef LA_WITH_XDP while sp, dp and their complex companions emit nothing.
              A loop over the kinds above the one its enclosing loop iterates over declares
              those kinds, so it carries their guard rather than the enclosing one's.

  tables      include/la_blas_interfaces.fypp and include/la_lapack_interfaces.fypp gain, for
              every qp and every complex-qp entry, the xdp entry of the same routine.
"""

import argparse
import ast
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import la_kindmap as K

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Loop sets whose iterations are one precision each.  The test templates alias two of them.
KIND_SETS = ("ALL_KINDS_TYPES", "REAL_KINDS_TYPES", "CMPL_KINDS_TYPES",
             "RC_KINDS_TYPES", "CMPLX_KINDS_TYPES", "LA_REAL_KINDS", "LA_CMPL_KINDS")
# The umbrella tables loop over interface entries; the precision is the leading letter.
ENTRY_SET = "entries"
ENTRY_KEY = "specific[0]"
# A loop over the kinds above the one its enclosing loop iterates over declares those kinds, so
# it carries their guard rather than the enclosing one's.
UP_SET = "LA_UPS["

FOR = re.compile(r"^(\s*)#:\s*for\s+(.+?)\s+in\s+(.+?)\s*$")
ENDFOR = re.compile(r"^\s*#:\s*endfor\s*$")
ENTRY_LINE = re.compile(r"^    &  (\[.*\]), &$")


def guard_key(line):
    """The expression whose LA_GUARD entry fences one iteration of this loop, or None."""
    m = FOR.match(line)
    if not m:
        return None
    targets, source = m.group(2), m.group(3)
    if source == ENTRY_SET or source.startswith(UP_SET) or source in KIND_SETS:
        return ENTRY_KEY if source == ENTRY_SET else targets.split(",")[0].strip()
    return None


def fence(indent, key, closing):
    directive = "#endif" if closing else "#ifdef ${LA_GUARD[%s]}$" % key
    return ["%s#:if LA_GUARD[%s]" % (indent, key), directive, "%s#:endif" % indent]


def guard_loops(text):
    """Wrap the body of every outermost kind loop in the cpp fence of its own precision.

    The fence sits inside the blank lines a loop body opens and closes with, so that the runs of
    blank lines the unguarded instances emit stay adjacent and fprettify collapses them as before.
    """
    lines = text.split("\n")
    out, depth, pending = [], 0, []
    i = 0
    while i < len(lines):
        line = lines[i]
        if pending and ENDFOR.match(line) and depth == pending[-1][0]:
            indent, key = pending.pop()[1:]
            closing, blanks = fence(indent, key, True), []
            while out and not out[-1].strip():
                blanks.insert(0, out.pop())
            if out[-3:] != closing:
                out += closing
            out += blanks
            out.append(line)
            depth -= 1
            i += 1
            continue
        if ENDFOR.match(line):
            depth -= 1
            out.append(line)
            i += 1
            continue
        key = guard_key(line)
        m = FOR.match(line)
        if m:
            depth += 1
        if key is None or (pending and not m.group(3).startswith(UP_SET)):
            out.append(line)
            i += 1
            continue
        opening = fence(m.group(1), key, False)
        out.append(line)
        i += 1
        while i < len(lines) and not lines[i].strip():
            out.append(lines[i])
            i += 1
        if lines[i:i + 3] != opening:
            out += opening
        pending.append((depth, m.group(1), key))
    return "\n".join(out)


def mirror_entries(entries):
    """The xdp entries of a generic: one per entry that names qp, by the same kind algebra."""
    extra = []
    known = {e[0] for e in entries}
    for specific, internal, stub in entries:
        mirrored = K.mirror_name(internal)
        if mirrored == internal:
            continue
        new = (mirrored, mirrored, stub)
        if new[0] not in known:
            extra.append(new)
    return extra


def guard_tables(text):
    """Add the xdp entry of every qp entry, keeping each list sorted by specific name."""
    out = []
    for line in text.split("\n"):
        m = ENTRY_LINE.match(line)
        if not m:
            out.append(line)
            continue
        entries = [tuple(e) for e in ast.literal_eval(m.group(1))]
        extra = mirror_entries(entries)
        if extra:
            entries = sorted(entries + extra, key=lambda e: e[0])
        out.append("    &  %r, &" % (entries,))
    return "\n".join(out)


def templates():
    for folder in ("src", "test"):
        base = os.path.join(ROOT, "fypp", folder)
        for name in sorted(os.listdir(base)):
            if name.endswith(".fypp"):
                yield os.path.join(base, name)


def tables():
    for name in ("la_blas_interfaces.fypp", "la_lapack_interfaces.fypp"):
        yield os.path.join(ROOT, "include", name)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true",
                        help="do not write; fail if a file is not already transformed")
    args = parser.parse_args()

    stale, written = [], 0
    for path, transform in ([(p, guard_loops) for p in templates()]
                            + [(p, guard_tables) for p in tables()]):
        with open(path) as fid:
            text = fid.read()
        new = transform(text)
        if new == text:
            continue
        if args.check:
            stale.append(os.path.relpath(path, ROOT))
            continue
        with open(path, "w") as fid:
            fid.write(new)
        written += 1

    if args.check:
        for name in stale:
            print("NOT GUARDED  " + name, file=sys.stderr)
        print("checked %d files, %d not guarded" % (len(list(templates())) + 2, len(stale)))
        return 1 if stale else 0
    print("guarded %d files" % written)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Kind-role mapping shared by templatize.py and check_generated.py.

A BLAS/LAPACK routine body is written for one precision.  `normalize` rewrites it into a
kind-neutral form by replacing every token that encodes a precision with a role placeholder:

    @RI@  own initial (s/d/q for a real routine, c/z/w for a complex one)
    @CI@  companion initial, the other type class at the same precision
    @RK@  own kind (sp/dp/qp)
    @RLI@ @RLK@ @CLI@   the same three, one precision down
    @RUI@ @RUK@ @CUI@   the same three, one precision up
    @PREC@ @PRECU@      the word a doc comment uses for the own precision, lower and upper case

Two bodies of the same routine at different precisions normalize to the same text unless they
really differ.  The placeholders map one to one onto the fypp loop variables of LA_REAL_KINDS
and LA_CMPL_KINDS in include/common.fypp.
"""

import collections
import os
import re

RI = ["s", "d", "q"]
CI = ["c", "z", "w"]
KND = ["sp", "dp", "qp"]
LETTERS = RI + CI
INDEX = {k: i for i, k in enumerate(RI)}
INDEX.update({k: i for i, k in enumerate(CI)})

PLACEHOLDERS = {
    "@RI@": "ri", "@CI@": "ci", "@RK@": "rk",
    "@RLI@": "ril", "@RLK@": "rkl", "@CLI@": "cil",
    "@RUI@": "riu", "@RUK@": "rku", "@CUI@": "ciu",
}

# Routine stems that carry no precision at all.  `iladiag` reads as `ila<d>iag` to the
# ila?lc/ila?lr rule below, which is what produced the removed la_ilaqiag copy.
KINDFREE = set("""lsame lsamen xerbla xerbla_array ilaenv ilaenv2stage ilatrans ilauplo
ilaprec ieeeck iparmq iparam2stage chla_transtype iladiag""".split())

PRECISION = {"sp": "single", "dp": "double", "qp": "quad"}

TWO_RC = ("asum", "nrm2", "sum1")      # scasum, dznrm2, ...:  real result over complex data
TWO_CR = ("rot", "scal", "rscl")       # csrot, zdscal, ...:   complex data, real scalar
MIXED = ("gesv", "posv", "dot")        # dsgesv, zcgesv, dsdot: two precisions


def roles(letter):
    """Placeholder map for a body whose own kind letter is `letter`.

    @R..I@ always stands for the initial of the routine's own type class and @C..I@ for the
    companion class, so a complex body maps c/z/w onto @R..I@ exactly as a real body maps s/d/q.
    """
    i = INDEX[letter]
    own, other = (CI, RI) if letter in CI else (RI, CI)
    label = {0: "", -1: "L", 1: "U", -2: "2L", 2: "2U"}
    out = {}
    for j in range(3):
        d = j - i
        if d not in label:
            continue
        out[own[j]] = "@R%sI@" % label[d]
        out[other[j]] = "@C%sI@" % label[d]
        out[KND[j]] = "@R%sK@" % label[d]
    return out


def map_name(stem, rl):
    """Normalize a routine stem, i.e. the text after the `la_` prefix. Case is preserved."""
    low = stem.lower()
    if low in KINDFREE:
        return stem
    sub = lambda ch: rl.get(ch, ch)
    for pattern, build in (
        (r"(lapack|blas)_([sdqczw])", lambda m: m.group(1) + "_" + sub(m.group(2))),
        (r"(selctg|select)_([sdqczw])", lambda m: m.group(1) + "_" + sub(m.group(2))),
        (r"i([sdqczw])(amax|max1)", lambda m: "i" + sub(m.group(1)) + m.group(2)),
        (r"ila([sdqczw])(lc|lr|iag)", lambda m: "ila" + sub(m.group(1)) + m.group(2)),
        (r"([sdq])([czw])(%s)" % "|".join(TWO_RC),
         lambda m: sub(m.group(1)) + sub(m.group(2)) + m.group(3)),
        (r"([czw])([sdq])(%s)" % "|".join(TWO_CR),
         lambda m: sub(m.group(1)) + sub(m.group(2)) + m.group(3)),
        (r"([sdqczw])(lag2|lat2)([sdqczw])",
         lambda m: sub(m.group(1)) + m.group(2) + sub(m.group(3))),
        (r"([sdq])([sdq])(%s)" % "|".join(MIXED),
         lambda m: sub(m.group(1)) + sub(m.group(2)) + m.group(3)),
        (r"([czw])([czw])(%s)" % "|".join(MIXED),
         lambda m: sub(m.group(1)) + sub(m.group(2)) + m.group(3)),
        (r"([sdqczw])(.*)", lambda m: sub(m.group(1)) + m.group(2)),
    ):
        m = re.fullmatch(pattern, low)
        if m:
            return build(m)
    return stem


def routine_names(paths):
    """Every `la_*` routine defined in the given Fortran sources, lower case, without `la_`."""
    pattern = re.compile(r"(?m)^\s*end\s+(?:subroutine|function)\s+la_([a-z0-9_]+)")
    names = set()
    for path in paths:
        with open(path, errors="replace") as fid:
            names.update(pattern.findall(fid.read()))
    return names


IMPORT = re.compile(r"(?m)^[ \t]*import\b.*$")


def normalize(text, letter, upper_bases=frozenset()):
    """Rewrite a body written for kind `letter` into placeholder form.

    An `import` statement of an interface body names the host entities the interface needs, not
    the precision it is written for, so it is carried through untouched.
    """
    kept = IMPORT.findall(text)
    if kept:
        text = IMPORT.sub("@IMPORT@", text)
    rl = roles(letter)
    text = re.sub(r"\bla_([A-Za-z0-9_]+)", lambda m: "la_" + map_name(m.group(1), rl), text)

    def upper(m):
        s = m.group(0)
        if s.lower() not in upper_bases:
            return s
        return map_name(s.lower(), rl).upper()

    text = re.sub(r"\b[A-Z][A-Z0-9_]{2,}\b", upper, text)
    for kind in KND:
        ph = rl.get(kind)
        if ph is None:
            continue
        text = re.sub(r"_%s\b" % kind, "_" + ph, text)
        text = re.sub(r"(?<![A-Za-z0-9_@])%s(?![A-Za-z0-9_])" % kind, ph, text)
    own = PRECISION[KND[INDEX[letter]]]
    text = re.sub(r"\b%s(?=[- ]precision\b)" % own, "@PREC@", text)
    text = re.sub(r"\b%s(?=[- ]PRECISION\b)" % own.upper(), "@PRECU@", text)
    if kept:
        lines = iter(kept)
        text = re.sub("@IMPORT@", lambda m: next(lines), text)
    return text


def strip_comment(line):
    """Drop a trailing `!` comment, ignoring `!` inside a character literal."""
    out, quoted = [], False
    for ch in line:
        if ch == "'":
            quoted = not quoted
        if ch == "!" and not quoted:
            break
        out.append(ch)
    return "".join(out)


WHITESPACE = re.compile(r"[ \t]+")


def ws_norm(text):
    """Collapse runs of blanks and drop empty lines, so indentation cannot mask a match."""
    out = []
    for line in text.split("\n"):
        line = WHITESPACE.sub(" ", line).strip()
        if line:
            out.append(line)
    return "\n".join(out)


TOKEN = re.compile(r"@?[A-Za-z_][A-Za-z0-9_]*@?|\d+\.?\d*(?:[eEdD][-+]?\d+)?|\S")


def tokens(text):
    return TOKEN.findall(text)


def letter_renames(donor_text, probe_text, pairs):
    """Identifiers that differ between two normalized bodies only by a leading kind letter.

    `pairs` maps (donor letter, probe letter) to the placeholder that stands for the pair, e.g.
    {("d", "s"): "@RI@", ("z", "c"): "@CI@"} for a real routine templated from the d body.
    Returns {donor identifier: placeholder + remainder}.
    """
    import difflib
    a, b = tokens(donor_text), tokens(probe_text)
    out = {}
    matcher = difflib.SequenceMatcher(a=a, b=b, autojunk=False)
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag != "replace" or (i2 - i1) != (j2 - j1):
            continue
        for x, y in zip(a[i1:i2], b[j1:j2]):
            if x == y or len(x) < 2 or len(y) < 2 or x[1:] != y[1:]:
                continue
            ph = pairs.get((x[0], y[0]))
            if ph is not None:
                out[x] = ph + x[1:]
    return out


def apply_renames(text, renames):
    if not renames:
        return text
    pattern = re.compile(r"\b(%s)\b" % "|".join(sorted(renames, key=len, reverse=True)))
    return pattern.sub(lambda m: renames[m.group(1)], text)


END_ROUTINE = re.compile(r"(?m)^\s*end\s+(subroutine|function)\s+la_([a-z0-9_]+)\s*$")


def split_routines(text):
    """Split a module body into `name -> (chunk, order)`, each chunk keeping its doc comment."""
    parts = text.split("\n     contains\n", 1)
    body = parts[1] if len(parts) > 1 else text
    out = collections.OrderedDict()
    pos, order = 0, 0
    for m in END_ROUTINE.finditer(body):
        out[m.group(2)] = (body[pos:m.end()], order)
        pos = m.end()
        order += 1
    return out


ABSTRACT = re.compile(r"(?ms)^([ \t]*abstract interface[ \t]*\n)(.*?)"
                      r"(^[ \t]*end interface[ \t]*$)")


def split_preamble(text):
    """Interface bodies declared before `contains`: `name -> (chunk, order)`.

    `abstract interface` blocks are part of a module's declaration section, so `split_routines`
    never sees them, yet their bodies carry a precision exactly as a routine does.
    """
    head = text.split("\n     contains\n", 1)[0]
    out = collections.OrderedDict()
    order = 0
    for block in ABSTRACT.finditer(head):
        body, pos = block.group(2), 0
        for m in END_ROUTINE.finditer(body):
            out[m.group(2)] = (body[pos:m.end()], order)
            pos = m.end()
            order += 1
    return out


def split_file(path):
    with open(path, errors="replace") as fid:
        return split_routines(fid.read())


def to_fypp(text):
    """Turn placeholder text into fypp: `@RI@gemm` -> `${ri}$gemm`, `@RI@GEMM` -> `${ri.upper()}$GEMM`."""
    text = text.replace("@PREC@", "${LA_PRECISION[rk]}$")
    text = text.replace("@PRECU@", "${LA_PRECISION[rk].upper()}$")
    run = re.compile(r"((?:@[A-Z0-9]+@)+)([A-Za-z0-9_]*)")

    def repl(m):
        marks = re.findall(r"@[A-Z0-9]+@", m.group(1))
        tail = m.group(2)
        upper = bool(re.search(r"[A-Z]", tail)) or (tail == "" and _upper_context(m))
        parts = []
        for mark in marks:
            var = PLACEHOLDERS[mark]
            parts.append("${%s.upper()}$" % var if upper else "${%s}$" % var)
        return "".join(parts) + tail

    def _upper_context(m):
        before = m.string[max(0, m.start() - 1):m.start()]
        return before.isupper()

    return run.sub(repl, text)

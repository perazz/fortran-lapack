#!/usr/bin/env python3
"""Kind-role mapping shared by templatize.py and check_generated.py.

A BLAS/LAPACK routine body is written for one precision.  `normalize` rewrites it into a
kind-neutral form by replacing every token that encodes a precision with a role placeholder:

    @RI@  own initial (s/d/x/q for a real routine, c/z/y/w for a complex one)
    @CI@  companion initial, the other type class at the same precision
    @RK@  own kind (sp/dp/xdp/qp)
    @RLI@ @RLK@ @CLI@   the same three, one precision down
    @RUI@ @RUK@ @CUI@   the same three, one precision up
    @PREC@ @PRECU@      the word a doc comment uses for the own precision, lower and upper case
    @PRECL@             the word a doc comment uses for the precision one below

Two bodies of the same routine at different precisions normalize to the same text unless they
really differ.  The placeholders map one to one onto the fypp loop variables of LA_REAL_KINDS
and LA_CMPL_KINDS in include/common.fypp.
"""

import collections
import os
import re

RI = ["s", "d", "x", "q"]
CI = ["c", "z", "y", "w"]
KND = ["sp", "dp", "xdp", "qp"]
LETTERS = RI + CI
CLASS = "[%s]" % "".join(LETTERS)
REAL_CLASS = "[%s]" % "".join(RI)
CMPL_CLASS = "[%s]" % "".join(CI)
INDEX = {k: i for i, k in enumerate(RI)}
INDEX.update({k: i for i, k in enumerate(CI)})
KIND_OF = {letter: KND[INDEX[letter]] for letter in LETTERS}

# The precision below each kind, and the precisions above it.  xdp is a side branch: it sits
# above dp and has nothing above it, so qp keeps dp below it and the mixed-precision routines of
# qp keep the names they have, while dp has two kinds above it.
DOWN = {"sp": None, "dp": "sp", "xdp": "dp", "qp": "dp"}
UP = {"sp": ["dp"], "dp": ["xdp", "qp"], "xdp": [], "qp": []}
# The xdp spelling of a qp routine, letter by letter.
MIRROR = {"q": "x", "w": "y"}

PLACEHOLDERS = {
    "@RI@": "ri", "@CI@": "ci", "@RK@": "rk",
    "@RLI@": "ril", "@RLK@": "rkl", "@CLI@": "cil",
    "@RUI@": "riu", "@RUK@": "rku", "@CUI@": "ciu",
}

# Routine stems that carry no precision at all.  `iladiag` reads as `ila<d>iag` to the
# ila?lc/ila?lr rule below, which is what produced the removed la_ilaqiag copy.
KINDFREE = set("""lsame lsamen xerbla xerbla_array ilaenv ilaenv2stage ilatrans ilauplo
ilaprec ieeeck iparmq iparam2stage chla_transtype iladiag""".split())

PRECISION = {"sp": "single", "dp": "double", "xdp": "extended", "qp": "quad"}

TWO_RC = ("asum", "nrm2", "sum1")      # scasum, dznrm2, ...:  real result over complex data
TWO_CR = ("rot", "scal", "rscl")       # csrot, zdscal, ...:   complex data, real scalar
MIXED = ("gesv", "posv", "dot")        # dsgesv, zcgesv, dsdot: two precisions


def roles(letter):
    """Placeholder map for a body whose own kind letter is `letter`.

    @R..I@ always stands for the initial of the routine's own type class and @C..I@ for the
    companion class, so a complex body maps c/z/y/w onto @R..I@ exactly as a real body maps
    s/d/x/q.  Only the kind itself and its two neighbours have a role; a precision further away
    keeps its own letter.
    """
    i = INDEX[letter]
    own, other = (CI, RI) if letter in CI else (RI, CI)
    out = {}
    below = [DOWN[KND[i]]] if DOWN[KND[i]] else []
    for label, kinds in (("", [KND[i]]), ("L", below), ("U", UP[KND[i]])):
        for kind in kinds:
            j = KND.index(kind)
            out[own[j]] = "@R%sI@" % label
            out[other[j]] = "@C%sI@" % label
            out[kind] = "@R%sK@" % label
    return out


MODULE_NAME = re.compile(r"(lapack|blas)_%s" % CLASS)


def map_name(stem, rl, known=None):
    r"""Normalize a routine stem, i.e. the text after the `la_` prefix. Case is preserved.

    `known` is the set of routine names of the tree.  Without it every `la_` token is treated as
    a routine, which turns the module name in `\see la_constants.f90` into a `c` routine.
    """
    low = stem.lower()
    if low in KINDFREE:
        return stem
    if known and low not in known and not MODULE_NAME.fullmatch(low):
        return stem
    sub = lambda ch: rl.get(ch, ch)
    for pattern, build in (
        (r"(lapack|blas)_(%s)" % CLASS, lambda m: m.group(1) + "_" + sub(m.group(2))),
        (r"(selctg|select)_(%s)" % CLASS, lambda m: m.group(1) + "_" + sub(m.group(2))),
        (r"i(%s)(amax|max1)" % CLASS, lambda m: "i" + sub(m.group(1)) + m.group(2)),
        (r"ila(%s)(lc|lr|iag)" % CLASS, lambda m: "ila" + sub(m.group(1)) + m.group(2)),
        (r"(%s)(%s)(%s)" % (REAL_CLASS, CMPL_CLASS, "|".join(TWO_RC)),
         lambda m: sub(m.group(1)) + sub(m.group(2)) + m.group(3)),
        (r"(%s)(%s)(%s)" % (CMPL_CLASS, REAL_CLASS, "|".join(TWO_CR)),
         lambda m: sub(m.group(1)) + sub(m.group(2)) + m.group(3)),
        (r"(%s)(lag2|lat2)(%s)" % (CLASS, CLASS),
         lambda m: sub(m.group(1)) + m.group(2) + sub(m.group(3))),
        (r"(%s)(%s)(%s)" % (REAL_CLASS, REAL_CLASS, "|".join(MIXED)),
         lambda m: sub(m.group(1)) + sub(m.group(2)) + m.group(3)),
        (r"(%s)(%s)(%s)" % (CMPL_CLASS, CMPL_CLASS, "|".join(MIXED)),
         lambda m: sub(m.group(1)) + sub(m.group(2)) + m.group(3)),
        (r"(%s)(.*)" % CLASS, lambda m: sub(m.group(1)) + m.group(2)),
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
    text = re.sub(r"\bla_([A-Za-z0-9_]+)",
                  lambda m: "la_" + map_name(m.group(1), rl, upper_bases), text)

    def upper(m):
        s = m.group(0)
        if s.lower() not in upper_bases:
            return s
        return map_name(s.lower(), rl, upper_bases).upper()

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


def lower_precision_word(text, letter):
    """Map the doc word of the precision one below `letter` onto @PRECL@.

    A routine that exists nowhere without a precision below it names that precision in its doc.
    Only the lower-case spaced form is templated: that is the one the per-kind copies rewrote,
    while the hyphenated form and the upper-case doc block keep the reference spelling.
    """
    below = DOWN[KND[INDEX[letter]]]
    if below is None:
        return text
    return re.sub(r"\b%s(?= precision\b)" % PRECISION[below], "@PRECL@", text)


def mirror_name(stem):
    """The xdp spelling of a qp routine stem: q becomes x and w becomes y, kind letters only."""
    return map_name(stem, MIRROR)


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


IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def letter_renames(donor_text, probe_text, pairs):
    """Identifiers that differ between two normalized bodies in one kind letter only.

    `pairs` maps (donor letter, probe letter) to the placeholder that stands for the pair, e.g.
    {("d", "s"): "@RI@", ("z", "c"): "@CI@"} for a real routine templated from the d body.  The
    letter is usually the first, as in `dnrm2` against `snrm2`, but the generator that wrote the
    q and w copies also rewrote it inside a local name, leaving `symb_wero` for `symb_zero`.
    Returns {donor identifier: placeholder in place of that letter}.
    """
    import difflib
    a, b = tokens(donor_text), tokens(probe_text)
    out = {}
    matcher = difflib.SequenceMatcher(a=a, b=b, autojunk=False)
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag != "replace" or (i2 - i1) != (j2 - j1):
            continue
        for x, y in zip(a[i1:i2], b[j1:j2]):
            if x == y or len(x) < 2 or len(x) != len(y) or not IDENTIFIER.fullmatch(x):
                continue
            at = [k for k in range(len(x)) if x[k] != y[k]]
            if len(at) != 1:
                continue
            ph = pairs.get((x[at[0]], y[at[0]]))
            if ph is not None:
                out[x] = x[:at[0]] + ph + x[at[0] + 1:]
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
    text = text.replace("@PRECL@", "${LA_PRECISION[rkl]}$")
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
        """A mark with no tail continues the identifier before it, if there is one."""
        head = re.search(r"[A-Za-z0-9_]+$", m.string[:m.start()])
        word = head.group(0) if head else ""
        return any(c.isupper() for c in word) and not any(c.islower() for c in word)

    return run.sub(repl, text)

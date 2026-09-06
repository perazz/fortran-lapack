#!/usr/bin/env python3
"""Turn the committed per-kind Fortran into kind-templated fypp topic modules.

    python3 scripts/templatize.py --library lapack --module la_lapack_solve_aux ...
                                                                         # write to the staging dir
    python3 scripts/templatize.py --library lapack --apply --module ...  # staging -> fypp/src
    python3 scripts/templatize.py --library lapack --extract --module ...# drop the converted
                                                                         # routines from the
                                                                         # per-kind sources
    python3 scripts/templatize.py --rename                               # renamed call sites
    python3 scripts/templatize.py --blas-interfaces                      # umbrella data table

--library selects the per-kind sources: src/la_blas_{s,d,q,c,z,w}.f90 plus src/la_blas_aux.f90,
or the same six names for LAPACK plus src/la_lapack_aux.f90.  The routine-to-module assignment
comes from scripts/la_modules.tsv.  Bodies are taken from the donor kinds (d for real routines,
z for complex ones) and rewritten into placeholder form by scripts/la_kindmap.py; the s and c
bodies of the same routine decide whether one template can serve every precision or whether the
routine needs an `#:if rk == "sp"` guard.  Bodies declared in an `abstract interface` block of a
module preamble template the same way and are emitted back into the preamble.

Conversion never writes into fypp/src directly.  It writes to the staging directory and --apply
copies from there, so that a guard minimised by hand in fypp/src is not silently overwritten by
a later run.  Re-run --apply on a module only when you mean to discard those edits.
"""

import argparse
import collections
import csv
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import la_kindmap as K

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAGING = os.path.join(ROOT, "build", "templatize")

REAL_KINDS = ["sp", "dp", "qp"]
REAL_INITIALS = ["s", "d", "q"]
CMPL_INITIALS = ["c", "z", "w"]

# One row of LA_REAL_KINDS / LA_CMPL_KINDS as include/common.fypp builds it.
def kind_rows(complex_class):
    own = CMPL_INITIALS if complex_class else REAL_INITIALS
    other = REAL_INITIALS if complex_class else CMPL_INITIALS
    rows = []
    for i, rk in enumerate(REAL_KINDS):
        down = i - 1 if i > 0 else None
        up = i + 1 if i < len(REAL_KINDS) - 1 else None
        rows.append({
            "rk": rk, "ri": own[i], "ci": other[i],
            "rkl": REAL_KINDS[down] if down is not None else None,
            "ril": own[down] if down is not None else None,
            "cil": other[down] if down is not None else None,
            "rku": REAL_KINDS[up] if up is not None else None,
            "riu": own[up] if up is not None else None,
            "ciu": other[up] if up is not None else None,
        })
    return rows


PLACEHOLDER_FIELD = {
    "@RI@": "ri", "@CI@": "ci", "@RK@": "rk",
    "@RLI@": "ril", "@RLK@": "rkl", "@CLI@": "cil",
    "@RUI@": "riu", "@RUK@": "rku", "@CUI@": "ciu",
}
DOWN_MARKS = ("@RLI@", "@RLK@", "@CLI@")
UP_MARKS = ("@RUI@", "@RUK@", "@CUI@")

COMPANION = {"s": "c", "d": "z", "q": "w", "c": "s", "z": "d", "w": "q"}
KIND_OF = {"s": "sp", "d": "dp", "q": "qp", "c": "sp", "z": "dp", "w": "qp"}


def _sources(prefix):
    return ["src/%s_%s.f90" % (prefix, k) for k in "sdqczw"] + ["src/%s_aux.f90" % prefix]


# One entry per library the converter reads: the per-kind sources, the donor body per type class
# (`d` and `z`, the most complete), the body divergence is probed against, the module that holds
# the helpers, and the aux routines whose name is not "<initial><stem>" and whose donor therefore
# cannot be built by prefixing the donor letter.
LIBRARIES = {
    "blas": {
        "prefix": "la_blas",
        "sources": _sources("la_blas"),
        "donor": {"real": ("d", "src/la_blas_d.f90"), "complex": ("z", "src/la_blas_z.f90")},
        "probe": {"real": ("s", "src/la_blas_s.f90"), "complex": ("c", "src/la_blas_c.f90")},
        "aux_module": "la_blas_aux",
        "aux_donor": {("iamax", "real"): "idamax", ("iamax", "complex"): "izamax",
                      ("cabs1", "real"): "dcabs1"},
    },
    "lapack": {
        "prefix": "la_lapack",
        "sources": _sources("la_lapack"),
        "donor": {"real": ("d", "src/la_lapack_d.f90"), "complex": ("z", "src/la_lapack_z.f90")},
        "probe": {"real": ("s", "src/la_lapack_s.f90"), "complex": ("c", "src/la_lapack_c.f90")},
        "aux_module": "la_lapack_aux",
        "aux_donor": {("ilalc", "real"): "iladlc", ("ilalc", "complex"): "ilazlc",
                      ("ilalr", "real"): "iladlr", ("ilalr", "complex"): "ilazlr",
                      ("imax1", "complex"): "izmax1",
                      ("select", "real"): "select_d", ("select", "complex"): "select_z",
                      ("selctg", "real"): "selctg_d", ("selctg", "complex"): "selctg_z"},
    },
}

MODULE_DOC = {
    "la_blas_aux": "BLAS helpers: character comparison, error reporting, index of maximum",
    "la_blas_level1": "BLAS level 1: vector operations",
    "la_blas_level2_ban": "BLAS level 2: banded matrix-vector operations",
    "la_blas_level2_gen": "BLAS level 2: general matrix-vector operations and rank updates",
    "la_blas_level2_pac": "BLAS level 2: packed and symmetric-banded matrix-vector operations",
    "la_blas_level2_sym": "BLAS level 2: symmetric matrix-vector operations",
    "la_blas_level2_tri": "BLAS level 2: triangular matrix-vector operations",
    "la_blas_level3_gen": "BLAS level 3: general and Hermitian matrix-matrix operations",
    "la_blas_level3_sym": "BLAS level 3: symmetric matrix-matrix operations",
    "la_blas_level3_tri": "BLAS level 3: triangular matrix-matrix operations",
    "la_lapack_aux": "LAPACK helpers: environment enquiry, character decoding, index scans",
    "la_lapack_auxiliary": "LAPACK auxiliary: machine parameters, safe division, band scaling",
    "la_lapack_blas_like_base": "BLAS-like base: copy, precision conversion, random and packed"
                                " storage",
    "la_lapack_blas_like_l1": "BLAS-like level 1: scaling, conjugation, sums of squares, sorting",
    "la_lapack_blas_like_l2": "BLAS-like level 2: matrix-vector products, scaling, rank updates",
    "la_lapack_blas_like_l3": "BLAS-like level 3: rank-k updates and solves in RFP storage",
    "la_lapack_blas_like_mnorm": "BLAS-like matrix norms",
    "la_lapack_blas_like_scalar": "BLAS-like scalar: complex division, Pythagorean sums, NaN"
                                  " tests",
    "la_lapack_givens_jacobi_rot": "Givens and Jacobi plane rotations",
    "la_lapack_householder_reflectors": "Householder reflectors: generation, blocking,"
                                        " application",
    "la_lapack_solve_aux": "Linear solve helpers: condition estimation, componentwise backward"
                           " error",
}


CONSTANTS = ("negone zero half one two three four eight ten czero chalf cone cnegone "
             "maxexp minexp rradix ulp eps safmin safmax smlnum bignum rtmin rtmax "
             "tsml tbig ssml sbig").split()
SIGNATURE = re.compile(r"^(\s*)\S.*\b(?:subroutine|function)\s+la_")
DECLARATION = re.compile(r"^\s*(?:real|complex|integer|logical|character)\s*(?:\([^)]*\))?"
                         r"\s*(?:,[^:]*)?::\s*(.*)$")
WORD = re.compile(r"[a-z_][a-z0-9_]*")


def _logical_lines(text):
    """Comment-free lines with `&` continuations folded into the line they continue."""
    out = []
    for line in text.split("\n"):
        line = K.strip_comment(line).rstrip()
        if out and out[-1].endswith("&"):
            out[-1] = out[-1][:-1] + line.strip()
        else:
            out.append(line)
    return out


def constants_used(text):
    """The constants of la_constants_<kind> a body reads, and the ones it declares itself."""
    lines = _logical_lines(text)
    declared = set()
    for line in lines:
        m = DECLARATION.match(line)
        if m:
            declared |= set(WORD.findall(m.group(1)))
    referenced = set(WORD.findall("\n".join(lines)))
    return ({c for c in CONSTANTS if c in referenced and c not in declared},
            {c for c in CONSTANTS if c in declared})


def needs_constants(record):
    """The import list for la_constants_<kind>: the names to bring in, and whether to name them.

    A LAPACK body routinely declares a local `safmin` or `eps` while reading `one` and `czero`
    from the module-level block, and a local declaration of a use-associated name is an error, so
    the import has to name what it brings in whenever the body shadows any of the constants.
    """
    wanted, shadowed = constants_used(record["text"])
    if record.get("sp_text"):
        other_wanted, other_shadowed = constants_used(record["sp_text"])
        wanted |= other_wanted
        shadowed |= other_shadowed
    clash = wanted & shadowed
    if clash:
        raise SystemExit("templatize: %s declares %s in one precision and reads it from "
                         "la_constants in the other" % (record["donor"], sorted(clash)))
    names = [c for c in CONSTANTS if c in wanted]
    return names, bool(shadowed)


def insert_use(text, kind, names):
    """Put `use la_constants_<kind>` right after the signature line."""
    lines = text.split("\n")
    for i, line in enumerate(lines):
        m = SIGNATURE.match(line)
        if not m:
            continue
        lines.insert(i + 1, m.group(1) + "   use la_constants_" + kind + names)
        return "\n".join(lines)
    raise SystemExit("templatize: no signature line in\n" + text[:200])


def read_modules(path):
    rows = []
    with open(path) as fid:
        for row in csv.DictReader(fid, delimiter="\t"):
            rows.append((row["stem"], row["class"], row["module"], int(row["level"])))
    return rows


def donor_name(cfg, stem, cls, module):
    if module == cfg["aux_module"] and (stem, cls) in cfg["aux_donor"]:
        return cfg["aux_donor"][(stem, cls)]
    if cls == "kindfree":
        return stem
    return cfg["donor"][cls][0] + stem


def concrete_name(norm_name, row):
    """`@RI@@CI@asum` + a kind row -> `dzasum`; None when the row lacks the neighbour."""
    try:
        return re.sub(r"@[A-Z0-9]+@", lambda m: _need(row, m.group(0)), norm_name)
    except KeyError:
        return None


def _need(row, mark):
    value = row[PLACEHOLDER_FIELD[mark]]
    if value is None:
        raise KeyError(mark)
    return value


def read_source(rel, ref):
    """A donor source as `ref` holds it, or from the working tree when the ref does not have it.

    The ref comes first so that a second run reads the same donors as the first: --extract takes
    the converted routines out of the per-kind sources in the tree, and a template must not be
    rebuilt from what is left behind.
    """
    run = subprocess.run(["git", "-C", ROOT, "show", "%s:%s" % (ref, rel)],
                         capture_output=True, text=True)
    if run.returncode == 0:
        return run.stdout
    path = os.path.join(ROOT, rel)
    if os.path.exists(path):
        with open(path, errors="replace") as fid:
            return fid.read()
    raise SystemExit("templatize: %s is neither in %s nor in the tree" % (rel, ref))


def ref_sources(ref):
    """Every generated `src/la_*` file in `ref`, as `relative path -> text`."""
    run = subprocess.run(["git", "-C", ROOT, "ls-tree", "-r", "--name-only", ref, "src/"],
                         capture_output=True, text=True, check=True)
    out = {}
    for rel in run.stdout.split("\n"):
        stem, ext = os.path.splitext(os.path.basename(rel))
        if ext in (".f90", ".F90") and stem.startswith("la_"):
            out[rel] = read_source(rel, ref)
    return out


FOREIGN_USE = re.compile(r"(?m)^\s*(use\s*(?:,\s*intrinsic\s*::)?\s*"
                         r"([a-z_][a-z0-9_]*)\s*(?:,\s*only\s*:\s*(.*?))?)\s*$")
END_ROUTINE = re.compile(r"(?m)^\s*end\s+(?:subroutine|function)\s+la_([a-z0-9_]+)")


class Library:
    """The committed per-kind sources of one library, split into routines and indexed by name.

    `chunks` holds the routine bodies of the `contains` section and `interfaces` the bodies
    declared in an `abstract interface` block of a module preamble; both carry a precision and
    template the same way, but only the first kind is a routine to the rest of the tooling.
    `preamble` keeps the text that surrounds those interface bodies, so the templated module can
    reproduce it, and `foreign_uses` the imports of modules outside the library.
    """

    def __init__(self, cfg, ref="origin/main", upper_from=None):
        self.cfg = cfg
        self.chunks, self.order, self.origin = {}, {}, {}
        self.interfaces, self.iface_order, self.preamble = {}, {}, {}
        self.foreign_uses = []
        texts = {rel: read_source(rel, ref) for rel in cfg["sources"]}
        for rel, text in sorted(texts.items()):
            for name, (chunk, order) in K.split_routines(text).items():
                self.chunks[name] = chunk
                self.order[name] = order
                self.origin[name] = rel
            for name, (chunk, order) in K.split_preamble(text).items():
                self.interfaces[name] = chunk
                self.iface_order[name] = order
                self.origin[name] = rel
            head = text.split("\n     contains\n", 1)[0]
            block = K.ABSTRACT.search(head)
            if block:
                self.preamble[rel] = (_after_publics(head[:block.end(1)]), block.group(3))
            for whole, module, only in FOREIGN_USE.findall(head):
                if not module.startswith("la_"):
                    self.foreign_uses.append((whole, [w for w in WORD.findall(only or "")]))
        self.upper_bases = set()
        for text in list(texts.values()) + list((upper_from or {}).values()):
            self.upper_bases.update(END_ROUTINE.findall(text))

    def body(self, name):
        return self.chunks.get(name, self.interfaces.get(name))

    def normalized(self, name, letter):
        return K.normalize(self.body(name), letter, self.upper_bases)


def _after_publics(head):
    """The declaration text that follows the last `public ::` line of a module preamble."""
    last = None
    for m in re.finditer(r"(?m)^\s*public\s*::.*$", head):
        last = m
    return head[last.end():] if last else head


def templatize_routine(lib, stem, cls, module):
    """Return a record describing one templated routine."""
    dname = donor_name(lib.cfg, stem, cls, module)
    if lib.body(dname) is None:
        raise SystemExit("templatize: no donor body for %s (%s) in %s" % (stem, cls, module))
    preamble = dname in lib.interfaces
    record = {"stem": stem, "class": cls, "module": module, "donor": dname,
              "order": lib.iface_order[dname] if preamble else lib.order[dname],
              "preamble": preamble, "guard": None, "renames": {}, "diverges": False}
    if cls == "kindfree":
        record["text"] = lib.body(dname).strip("\n")
        record["norm_name"] = dname
        return _with_guard(record)

    dletter = _own_letter(dname)
    donor = lib.normalized(dname, dletter)
    norm_name = K.map_name(dname, K.roles(dletter))
    record["norm_name"] = norm_name

    # The single-precision body of the same routine decides whether one template serves all
    # precisions.  It may not exist (a mixed-precision routine has no sp shape).
    sp_row = kind_rows(cls == "complex")[0]
    pname = concrete_name(norm_name, sp_row)
    if pname is None or lib.body(pname) is None or pname == dname:
        record["text"] = donor.strip("\n")
        return _with_guard(record)

    pletter = _own_letter(pname)
    probe = lib.normalized(pname, pletter)
    renames = K.letter_renames(donor, probe,
                               {(dletter, pletter): "@RI@",
                                (COMPANION[dletter], COMPANION[pletter]): "@CI@"})
    donor_t = K.apply_renames(donor, renames)
    probe_t = K.apply_renames(probe, K.letter_renames(
        probe, donor, {(pletter, dletter): "@RI@",
                       (COMPANION[pletter], COMPANION[dletter]): "@CI@"}))

    record["renames"] = renames
    record["probe"] = pname
    if K.ws_norm(donor_t) == K.ws_norm(probe_t):
        record["text"] = donor_t.strip("\n")
    else:
        record["diverges"] = True
        record["text"] = donor_t.strip("\n")
        record["sp_text"] = probe_t.strip("\n")
    return _with_guard(record)


def _with_guard(record):
    marks = set(re.findall(r"@[A-Z0-9]+@", record["text"]))
    if record.get("sp_text"):
        marks |= set(re.findall(r"@[A-Z0-9]+@", record["sp_text"]))
    if marks & set(DOWN_MARKS):
        record["guard"] = "rkl is not None"
    elif marks & set(UP_MARKS):
        record["guard"] = "rku is not None"
    names, shadowed = needs_constants(record)
    if names:
        kind = "${rk}$" if record["class"] != "kindfree" else KIND_OF[_own_letter(record["donor"])]
        only = ",only:" + ",".join(names) if shadowed else ""
        record["text"] = insert_use(record["text"], kind, only)
        if record.get("sp_text"):
            record["sp_text"] = insert_use(record["sp_text"], kind, only)
    return record


OWN_LETTER = (r"(?:selctg|select)_([sdqczw])", r"i([sdqczw])(?:amax|max1)",
              r"ila([sdqczw])(?:lc|lr|iag)", r"([sdqczw]).*")


def _own_letter(name):
    for pattern in OWN_LETTER:
        m = re.fullmatch(pattern, name)
        if m:
            return m.group(1)
    raise SystemExit("templatize: cannot tell the kind letter of %s" % name)


def minimize(sp_text, text):
    """Guard only the lines where the single-precision body differs from the donor body."""
    import difflib
    sp_lines, lines = sp_text.split("\n"), text.split("\n")
    matcher = difflib.SequenceMatcher(a=[l.strip() for l in sp_lines],
                                      b=[l.strip() for l in lines], autojunk=False)
    out, hunks, guarded = [], 0, 0
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            out += lines[j1:j2]
            continue
        hunks += 1
        guarded += (i2 - i1) + (j2 - j1)
        if i2 > i1 and j2 > j1:
            out.append('#:if rk == "sp"')
            out += sp_lines[i1:i2]
            out.append("#:else")
            out += lines[j1:j2]
            out.append("#:endif")
        elif i2 > i1:
            out.append('#:if rk == "sp"')
            out += sp_lines[i1:i2]
            out.append("#:endif")
        else:
            out.append('#:if rk != "sp"')
            out += lines[j1:j2]
            out.append("#:endif")
    return "\n".join(out), hunks, guarded, len(lines)


def emit_routine(record):
    if record["class"] == "kindfree":
        return K.to_fypp(record["text"])
    loop = "LA_CMPL_KINDS" if record["class"] == "complex" else "LA_REAL_KINDS"
    lines = ["#:for rk, rt, ri, ci, rkl, ril, cil, rku, riu, ciu in " + loop]
    if record["guard"]:
        lines.append("#:if " + record["guard"])
    if record["diverges"]:
        body, hunks, guarded, total = minimize(record["sp_text"], record["text"])
        record["hunks"] = hunks
        record["guarded_lines"] = guarded
        record["total_lines"] = total
        lines.append(K.to_fypp(body))
    else:
        lines.append(K.to_fypp(record["text"]))
    if record["guard"]:
        lines.append("#:endif")
    lines.append("#:endfor")
    return "\n".join(lines)


def emit_public(records):
    """`public ::` lines, grouped by the loop that produces the names."""
    out = []
    for loop, cls in (("LA_REAL_KINDS", "real"), ("LA_CMPL_KINDS", "complex")):
        group = [r for r in records if r["class"] == cls]
        if not group:
            continue
        out.append("#:for rk, rt, ri, ci, rkl, ril, cil, rku, riu, ciu in " + loop)
        for r in group:
            if r["guard"]:
                out.append("#:if " + r["guard"])
            out.append("     public :: la_" + K.to_fypp(r["norm_name"]))
            if r["guard"]:
                out.append("#:endif")
        out.append("#:endfor")
    for r in records:
        if r["class"] == "kindfree":
            out.append("     public :: la_" + r["norm_name"])
    return out


CALL_RE = re.compile(r"\bla_([a-z0-9_]+)")


def module_uses(lib, records, name_to_module, module):
    """Modules the bodies reference, comments stripped before extracting the references.

    References are taken from the concrete donor sources, not from the templated text, where a
    callee name is broken up by placeholders.  A module outside the library is carried over from
    the header of the source the bodies came from, restricted to the ones they really name.
    """
    used, words = set(), set()
    for r in records:
        raw = lib.body(r["donor"]) + ("\n" + lib.body(r["probe"]) if r.get("probe") else "")
        text = "\n".join(K.strip_comment(line) for line in raw.split("\n"))
        words |= set(WORD.findall(text.lower()))
        for ref in CALL_RE.findall(text):
            owner = name_to_module.get(ref)
            if owner and owner != module:
                used.add(owner)
    foreign = [line for line, only in lib.foreign_uses if not only or (set(only) & words)]
    return ["use " + m for m in sorted(used)] + foreign


def emitted_names(cfg, stem, cls, module):
    """The concrete routine names the template of one row will emit, in kind order."""
    dname = donor_name(cfg, stem, cls, module)
    if cls == "kindfree":
        return [dname]
    norm = K.map_name(dname, K.roles(_own_letter(dname)))
    return [n for n in (concrete_name(norm, row) for row in kind_rows(cls == "complex")) if n]


def topic_owners(tree, modules):
    """`routine name -> module` for the topic modules the baseline tree already holds."""
    out = {}
    for rel, text in tree.items():
        stem = os.path.splitext(os.path.basename(rel))[0]
        if stem not in modules:
            continue
        for name in list(K.split_routines(text)) + list(K.split_preamble(text)):
            out[name] = stem
    return out


def build(args):
    cfg = LIBRARIES[args.library]
    rows = read_modules(args.modules)
    tree = ref_sources(args.baseline)
    lib = Library(cfg, args.baseline, upper_from=tree)
    mine = [r for r in rows if r[2].startswith(cfg["prefix"])]
    wanted = args.module or sorted({m for _, _, m, _ in mine})

    # Where every routine ends up: the topic modules the baseline already holds, then the ones
    # this library's rows describe, so that a `use` names the module a callee really lives in.
    name_to_module = topic_owners(tree, {m for _, _, m, _ in rows})
    emitted = {}
    for stem, cls, module, _ in mine:
        for name in emitted_names(cfg, stem, cls, module):
            name_to_module[name] = module
            emitted[name] = module

    report = []
    os.makedirs(args.stage, exist_ok=True)
    for module in wanted:
        members = [(s, c) for s, c, m, _ in rows if m == module]
        records = [templatize_routine(lib, stem, cls, module) for stem, cls in members]
        records.sort(key=lambda r: ({"real": 0, "complex": 1, "kindfree": 2}[r["class"]],
                                    r["order"]))
        uses = ["use la_constants"] + module_uses(lib, records, name_to_module, module)
        body = []
        body.append('#:include "common.fypp"')
        body.append("!> " + MODULE_DOC.get(module, module))
        body.append("module " + module)
        for use in uses:
            body.append("     " + use)
        body.append("     implicit none(type,external)")
        body.append("     private")
        body.append("")
        body.append("     public :: sp,dp,qp,lk,ilp")
        body += emit_public(records)
        body += emit_preamble(lib, [r for r in records if r["preamble"]])
        body.append("")
        body.append("     contains")
        for record in records:
            if record["preamble"]:
                continue
            body.append("")
            body.append(emit_routine(record))
        body.append("")
        body.append("end module " + module)
        body.append("")
        path = os.path.join(args.stage, module + ".fypp")
        with open(path, "w") as fid:
            fid.write("\n".join(body))
        for record in records:
            if record["diverges"]:
                report.append("%s\t%s\t%s\t%d guarded hunk(s), %d guarded lines of %d"
                              % (module, record["stem"], record["class"], record.get("hunks", 0),
                                 record.get("guarded_lines", 0), record.get("total_lines", 0)))
        print("staged %-36s %3d routines, %s"
              % (module, len(records), "; ".join(uses)))

    # Names that change: the committed source has one spelling, the template emits another.
    known = set(lib.chunks) | set(lib.interfaces)
    gone = sorted(known - set(emitted))
    fresh = sorted(set(emitted) - known)
    if gone or fresh:
        report.append("baseline names no longer emitted: " + ",".join(gone))
        report.append("names the templates add:          " + ",".join(fresh))
    if args.report:
        with open(args.report, "w") as fid:
            fid.write("\n".join(report) + "\n")
    for line in report:
        print("report: " + line)


def emit_preamble(lib, records):
    """The declaration-section text that surrounds the templated `abstract interface` bodies."""
    if not records:
        return []
    rel = lib.origin[records[0]["donor"]]
    head, foot = lib.preamble[rel]
    out = [""] + head.strip("\n").split("\n")
    for record in records:
        out.append(emit_routine(record))
    out.append(foot)
    return out


def apply_staged(args):
    for module in args.module:
        src = os.path.join(args.stage, module + ".fypp")
        dst = os.path.join(ROOT, "fypp", "src", module + ".fypp")
        shutil.copyfile(src, dst)
        print("applied %s -> %s" % (src, os.path.relpath(dst, ROOT)))


WITH_QP_OPEN = re.compile(r"(?m)^#!if WITH_QP[ \t]*$")
WITH_QP_CLOSE = re.compile(r"\A[ \t]*\n?#!endif[ \t]*\n")


def extract(args):
    """Delete the converted routines from the per-kind sources they came from.

    A routine the templates rename is deleted under the name the per-kind source gives it, and
    the inert `#!if WITH_QP` comment the quad copies are wrapped in goes with the body.
    """
    cfg = LIBRARIES[args.library]
    rows = read_modules(args.modules)
    lib = Library(cfg, args.baseline)
    renames = read_renames(os.path.join(ROOT, "scripts", "la_renames.tsv"))
    was = {new[3:]: old[3:] for old, new in renames.items()}
    targets = collections.defaultdict(list)
    for stem, cls, module, _ in rows:
        if args.module and module not in args.module:
            continue
        if not module.startswith(cfg["prefix"]):
            continue
        for name in emitted_names(cfg, stem, cls, module):
            name = name if name in lib.origin else was.get(name)
            rel = lib.origin.get(name)
            # A module that keeps its own name is replaced by its template, not extracted from.
            if rel and os.path.splitext(os.path.basename(rel))[0] != module:
                targets[rel].append(name)
    for rel, names in sorted(targets.items()):
        for path in (os.path.join(args.tree, rel),
                     os.path.join(args.tree, "fypp", rel.replace(".f90", ".fypp"))):
            if not os.path.exists(path):
                continue
            with open(path) as fid:
                text = fid.read()
            for name in names:
                chunks = K.split_routines(text)
                if name in chunks:
                    text = _drop_chunk(text, chunks[name][0])
                text = _drop_public(text, name)
            with open(path, "w") as fid:
                fid.write(text)
            print("extracted %d routines from %s" % (len(names), os.path.relpath(path, args.tree)))


def _drop_chunk(text, chunk):
    start = text.index(chunk)
    end = start + len(chunk)
    close = WITH_QP_CLOSE.match(text[end:])
    if close and WITH_QP_OPEN.search(chunk):
        end += close.end()
    return text[:start] + text[end:]


def _drop_public(text, name):
    return re.sub(r"(?m)^(?:#!if WITH_QP[ \t]*\n)? *public :: la_%s *\n(?:#!endif[ \t]*\n)?"
                  % re.escape(name), "", text)


def rename(args):
    """Apply scripts/la_renames.tsv to the sources that still spell the old specific names."""
    renames = read_renames(os.path.join(ROOT, "scripts", "la_renames.tsv"))
    if not renames:
        return
    pattern = re.compile(r"\b(%s)\b" % "|".join(sorted(renames, key=len, reverse=True)))
    for folder in ("src", os.path.join("fypp", "src")):
        base = os.path.join(args.tree, folder)
        for name in sorted(os.listdir(base)):
            path = os.path.join(base, name)
            with open(path) as fid:
                text = fid.read()
            new = pattern.sub(lambda m: renames[m.group(1)], text)
            if new == text:
                continue
            with open(path, "w") as fid:
                fid.write(new)
            print("renamed %d call sites in %s"
                  % (len(pattern.findall(text)), os.path.relpath(path, args.tree)))


IMPORT_LINE = "                    import sp,dp,qp,ilp,lk"
IMPORT_MASK = "                    import @IMPORTS@"


def read_renames(path):
    """`old_name -> new_name` for the specific names the templates spell differently."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fid:
        for row in csv.DictReader(fid, delimiter="\t"):
            if row["new_name"] != "REMOVED":
                out[row["old_name"]] = row["new_name"]
    return out


def blas_interfaces(args):
    """Build include/la_blas_interfaces.fypp from the committed umbrella."""
    lib = Library(LIBRARIES["blas"], args.baseline)
    renames = read_renames(os.path.join(ROOT, "scripts", "la_renames.tsv"))
    lines = read_source("src/la_blas.F90", args.baseline).split("\n")
    table, i = [], 0
    while i < len(lines):
        m = re.match(r"^          interface (\w+)$", lines[i])
        if not m:
            i += 1
            continue
        generic = m.group(1)
        j, doc = i - 1, []
        while j >= 0 and lines[j].startswith("          !>"):
            doc.insert(0, lines[j][len("          !> "):] if lines[j] != "          !>" else "")
            j -= 1
        k, entries, stubs = i + 1, [], {}
        while lines[k] != "          end interface " + generic:
            if lines[k] == "#ifdef LA_EXTERNAL_BLAS":
                j2, stub = k + 1, []
                while lines[j2] != "#else":
                    stub.append(lines[j2])
                    j2 += 1
                specific = re.match(r"\s*module procedure la_(\w+)$", lines[j2 + 1]).group(1)
                entries.append((specific, renames.get("la_" + specific, "la_" + specific)[3:],
                                True))
                cls = "real" if specific[0] in "sdq" else "complex"
                stubs.setdefault(cls, (specific, "\n".join(stub)))
                k = j2 + 3
            else:
                specific = re.match(r"\s*module procedure la_(\w+)$", lines[k]).group(1)
                entries.append((specific, renames.get("la_" + specific, "la_" + specific)[3:],
                                False))
                k += 1
        table.append((generic, doc, entries, stubs, lib))
        i = k + 1

    out = ["#:mute", "",
           "#! Generic BLAS interfaces of module la_blas, one entry per generic: name, the lines",
           "#! of its doc comment, (external specific name, module procedure, external stub) in",
           "#! file order, and one external-stub template per type class, with {ri}/{ci}/{rt}/{rk}",
           "#! and the neighbour fields of LA_BY_INITIAL as substitution fields.  Regenerate with",
           "#! `python3 scripts/templatize.py --blas-interfaces`.", "",
           "#:set LA_BLAS_INTERFACES = [ &"]
    for generic, doc, entries, stubs, _ in table:
        templates = {}
        for cls, (specific, text) in stubs.items():
            letter = specific[0]
            norm = K.normalize(text.replace(IMPORT_LINE, IMPORT_MASK), letter, lib.upper_bases)
            probe_letter = "s" if cls == "real" else "c"
            other = [(s, e) for s, _i, e in _same_class(entries, cls) if s[0] != letter]
            if other:
                pletter = other[0][0][0]
                pnorm = K.normalize(_stub_of(lines, generic, other[0][0]).replace(
                    IMPORT_LINE, IMPORT_MASK), pletter, lib.upper_bases)
                norm = K.apply_renames(norm, K.letter_renames(
                    norm, pnorm, {(letter, pletter): "@RI@",
                                  (COMPANION[letter], COMPANION[pletter]): "@CI@"}))
            templates[cls] = _to_format(norm.replace(IMPORT_MASK, IMPORT_LINE), cls)
        out.append("    & (%r, %r, &" % (generic, doc))
        out.append("    &  %r, &" % ([(s, i, bool(e)) for s, i, e in entries],))
        out.append("    &  {%s}), &"
                   % ", ".join("%r: %r" % (c, x) for c, x in sorted(templates.items())))
    out += ["    & ]", "", "#:endmute", ""]
    path = os.path.join(ROOT, "include", "la_blas_interfaces.fypp")
    with open(path, "w") as fid:
        fid.write("\n".join(out))
    print("wrote %s (%d generics)" % (os.path.relpath(path, ROOT), len(table)))


def _same_class(entries, cls):
    return [(s, i, e) for s, i, e in entries
            if e and (("real" if s[0] in "sdq" else "complex") == cls)]


def _stub_of(lines, generic, specific):
    i = lines.index("          interface " + generic)
    while lines[i] != "          end interface " + generic:
        if lines[i] == "#ifdef LA_EXTERNAL_BLAS":
            j, stub = i + 1, []
            while lines[j] != "#else":
                stub.append(lines[j])
                j += 1
            if lines[j + 1].strip() == "module procedure la_" + specific:
                return "\n".join(stub)
            i = j + 3
        else:
            i += 1
    raise KeyError(specific)


def _to_format(text, cls):
    """Placeholder text -> a python format string with {ri}/{ci}/{rk}/{rt} fields.

    Only the type of the routine's own class becomes {rt}; the companion type keeps its keyword
    and takes the kind from {rk}.
    """
    text = text.replace("{", "{{").replace("}", "}}")
    own = "complex" if cls == "complex" else "real"
    text = re.sub(r"\b%s\(@RK@\)" % own, "{rt}", text)
    return re.sub(r"@[A-Z0-9]+@", lambda m: "{%s}" % PLACEHOLDER_FIELD[m.group(0)], text)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--library", choices=sorted(LIBRARIES), default="blas",
                        help="which per-kind sources the donors come from")
    parser.add_argument("--modules", default=os.path.join(ROOT, "scripts", "la_modules.tsv"))
    parser.add_argument("--baseline", default="origin/main",
                        help="git ref the donor bodies come from once they leave the tree")
    parser.add_argument("--module", action="append", default=[])
    parser.add_argument("--stage", default=STAGING)
    parser.add_argument("--report", default=None)
    parser.add_argument("--apply", action="store_true", help="copy staged templates to fypp/src")
    parser.add_argument("--extract", action="store_true",
                        help="delete the converted routines from the per-kind sources")
    parser.add_argument("--tree", default=ROOT, help="tree --extract edits (a copy, for testing)")
    parser.add_argument("--blas-interfaces", action="store_true",
                        help="regenerate include/la_blas_interfaces.fypp")
    parser.add_argument("--rename", action="store_true",
                        help="apply scripts/la_renames.tsv to the call sites left behind")
    args = parser.parse_args()
    if args.blas_interfaces:
        return blas_interfaces(args)
    if args.apply:
        if not args.module:
            raise SystemExit("--apply needs at least one --module")
        return apply_staged(args)
    if args.extract:
        return extract(args)
    if args.rename:
        return rename(args)
    return build(args)


if __name__ == "__main__":
    sys.exit(main() or 0)

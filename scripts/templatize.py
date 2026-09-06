#!/usr/bin/env python3
"""Turn the committed per-kind Fortran into kind-templated fypp topic modules.

    python3 scripts/templatize.py --module la_blas_level1 --module ...   # write to the staging dir
    python3 scripts/templatize.py --apply --module la_blas_level1        # staging -> fypp/src
    python3 scripts/templatize.py --extract --module la_blas_level1      # drop the converted
                                                                         # routines from the
                                                                         # per-kind sources
    python3 scripts/templatize.py --blas-interfaces                      # umbrella data table

The routine-to-module assignment comes from scripts/la_modules.tsv.  Bodies are taken from the
donor kinds (d for real routines, z for complex ones) and rewritten into placeholder form by
scripts/la_kindmap.py; the s and c bodies of the same routine decide whether one template can
serve every precision or whether the routine needs an `#:if rk == "sp"` guard.

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

BLAS_SOURCES = ["src/la_blas_%s.f90" % k for k in "sdqczw"] + ["src/la_blas_aux.f90"]

# Donor bodies per type class, and the body used to detect divergence between precisions.
DONOR = {"real": ("d", "src/la_blas_d.f90"), "complex": ("z", "src/la_blas_z.f90")}
PROBE = {"real": ("s", "src/la_blas_s.f90"), "complex": ("c", "src/la_blas_c.f90")}
COMPANION = {"s": "c", "d": "z", "q": "w", "c": "s", "z": "d", "w": "q"}
KIND_OF = {"s": "sp", "d": "dp", "q": "qp", "c": "sp", "z": "dp", "w": "qp"}

# Routines of la_blas_aux whose name does not follow "<initial><stem>".
AUX_DONOR = {("iamax", "real"): "idamax", ("iamax", "complex"): "izamax",
             ("cabs1", "real"): "dcabs1"}

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


def needs_constants(record):
    """Whether the routine reads a constant of la_constants_<kind> that it does not declare."""
    lines = _logical_lines(record["text"] + "\n" + record.get("sp_text", ""))
    declared = set()
    for line in lines:
        m = DECLARATION.match(line)
        if m:
            declared |= set(WORD.findall(m.group(1)))
    referenced = set(WORD.findall("\n".join(lines)))
    wanted = [c for c in CONSTANTS if c in referenced and c not in declared]
    shadowed = [c for c in CONSTANTS if c in declared]
    if wanted and shadowed:
        raise SystemExit("templatize: %s declares %s locally but also needs %s from "
                         "la_constants" % (record["donor"], shadowed, wanted))
    return bool(wanted)


def insert_use(text, kind):
    """Put a bare `use la_constants_<kind>` right after the signature line."""
    lines = text.split("\n")
    for i, line in enumerate(lines):
        m = SIGNATURE.match(line)
        if not m:
            continue
        lines.insert(i + 1, m.group(1) + "   use la_constants_" + kind)
        return "\n".join(lines)
    raise SystemExit("templatize: no signature line in\n" + text[:200])


def read_modules(path):
    rows = []
    with open(path) as fid:
        for row in csv.DictReader(fid, delimiter="\t"):
            rows.append((row["stem"], row["class"], row["module"], int(row["level"])))
    return rows


def donor_name(stem, cls, module):
    if module == "la_blas_aux" and (stem, cls) in AUX_DONOR:
        return AUX_DONOR[(stem, cls)]
    if cls == "kindfree":
        return stem
    return DONOR[cls][0] + stem


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
    """A donor source from the working tree, or from `ref` once the conversion removed it."""
    path = os.path.join(ROOT, rel)
    if os.path.exists(path):
        with open(path, errors="replace") as fid:
            return fid.read()
    run = subprocess.run(["git", "-C", ROOT, "show", "%s:%s" % (ref, rel)],
                         capture_output=True, text=True)
    if run.returncode != 0:
        raise SystemExit("templatize: %s is neither in the tree nor in %s" % (rel, ref))
    return run.stdout


class Library:
    """The committed per-kind sources, split into routines and indexed by name."""

    def __init__(self, sources, ref="origin/main"):
        self.chunks = {}
        self.order = {}
        self.origin = {}
        texts = {rel: read_source(rel, ref) for rel in sources}
        for rel, text in texts.items():
            for name, (chunk, order) in K.split_routines(text).items():
                self.chunks[name] = chunk
                self.order[name] = order
                self.origin[name] = rel
        pattern = re.compile(r"(?m)^\s*end\s+(?:subroutine|function)\s+la_([a-z0-9_]+)")
        self.upper_bases = set()
        for text in texts.values():
            self.upper_bases.update(pattern.findall(text))

    def normalized(self, name, letter):
        return K.normalize(self.chunks[name], letter, self.upper_bases)


def templatize_routine(lib, stem, cls, module):
    """Return a record describing one templated routine."""
    dname = donor_name(stem, cls, module)
    if dname not in lib.chunks:
        raise SystemExit("templatize: no donor body for %s (%s) in %s" % (stem, cls, module))
    record = {"stem": stem, "class": cls, "module": module, "donor": dname,
              "order": lib.order[dname], "guard": None, "renames": {}, "diverges": False}
    if cls == "kindfree":
        record["text"] = lib.chunks[dname].strip("\n")
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
    if pname is None or pname not in lib.chunks or pname == dname:
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
    if needs_constants(record):
        kind = "${rk}$" if record["class"] != "kindfree" else KIND_OF[_own_letter(record["donor"])]
        record["text"] = insert_use(record["text"], kind)
        if record.get("sp_text"):
            record["sp_text"] = insert_use(record["sp_text"], kind)
    return record


OWN_LETTER = (r"i([sdqczw])(?:amax|max1)", r"ila([sdqczw])(?:lc|lr|iag)", r"([sdqczw]).*")


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
    """Modules referenced by the bodies, comments stripped before extracting the references.

    References are taken from the concrete donor sources, not from the templated text, where a
    callee name is broken up by placeholders.
    """
    used = set()
    for r in records:
        raw = lib.chunks[r["donor"]] + ("\n" + lib.chunks[r["probe"]] if r.get("probe") else "")
        text = "\n".join(K.strip_comment(line) for line in raw.split("\n"))
        for ref in CALL_RE.findall(text):
            owner = name_to_module.get(ref)
            if owner and owner != module:
                used.add(owner)
    return sorted(used)


def build(args):
    rows = read_modules(args.modules)
    lib = Library(BLAS_SOURCES, args.baseline)
    wanted = args.module or sorted({m for _, _, m, _ in rows if m.startswith("la_blas")})

    # Every concrete routine name the templates will emit, and where it ends up.
    name_to_module, emitted, renamed = {}, {}, []
    for stem, cls, module, _ in rows:
        if not module.startswith("la_blas"):
            continue
        record = templatize_routine(lib, stem, cls, module)
        if cls == "kindfree":
            name_to_module[record["norm_name"]] = module
            emitted[record["norm_name"]] = module
            continue
        for row in kind_rows(cls == "complex"):
            name = concrete_name(record["norm_name"], row)
            if name is None:
                continue
            name_to_module[name] = module
            emitted[name] = module

    report = []
    os.makedirs(args.stage, exist_ok=True)
    for module in wanted:
        members = [(s, c) for s, c, m, _ in rows if m == module]
        records = []
        for stem, cls in members:
            records.append(templatize_routine(lib, stem, cls, module))
        records.sort(key=lambda r: ({"real": 0, "complex": 1, "kindfree": 2}[r["class"]],
                                    r["order"]))
        uses = ["la_constants"] + module_uses(lib, records, name_to_module, module)
        body = []
        body.append('#:include "common.fypp"')
        body.append("!> " + MODULE_DOC.get(module, module))
        body.append("module " + module)
        for use in uses:
            body.append("     use " + use)
        body.append("     implicit none(type,external)")
        body.append("     private")
        body.append("")
        body.append("     public :: sp,dp,qp,lk,ilp")
        body += emit_public(records)
        body.append("")
        body.append("     contains")
        for record in records:
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
        print("staged %-24s %2d routines, uses %s" % (module, len(records), ",".join(uses)))

    # Names that change: the committed source has one spelling, the template emits another.
    gone = sorted(set(lib.chunks) - set(emitted))
    fresh = sorted(set(emitted) - set(lib.chunks))
    if gone or fresh:
        report.append("baseline names no longer emitted: " + ",".join(gone))
        report.append("names the templates add:          " + ",".join(fresh))
    if args.report:
        with open(args.report, "w") as fid:
            fid.write("\n".join(report) + "\n")
    for line in report:
        print("report: " + line)


def apply_staged(args):
    for module in args.module:
        src = os.path.join(args.stage, module + ".fypp")
        dst = os.path.join(ROOT, "fypp", "src", module + ".fypp")
        shutil.copyfile(src, dst)
        print("applied %s -> %s" % (src, os.path.relpath(dst, ROOT)))


def extract(args):
    """Delete the converted routines from the per-kind sources they came from."""
    rows = read_modules(args.modules)
    lib = Library(BLAS_SOURCES, args.baseline)
    targets = collections.defaultdict(list)
    for stem, cls, module, _ in rows:
        if args.module and module not in args.module:
            continue
        if not module.startswith("la_blas"):
            continue
        record = templatize_routine(lib, stem, cls, module)
        if cls == "kindfree":
            names = [record["norm_name"]]
        else:
            names = [concrete_name(record["norm_name"], row)
                     for row in kind_rows(cls == "complex")]
        for name in names:
            if name and name in lib.origin:
                targets[lib.origin[name]].append(name)
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
                    text = text.replace(chunks[name][0], "", 1)
                text = re.sub(r"(?m)^ *public :: la_%s *\n" % re.escape(name), "", text)
            with open(path, "w") as fid:
                fid.write(text)
            print("extracted %d routines from %s" % (len(names), os.path.relpath(path, args.tree)))


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
    lib = Library(BLAS_SOURCES, args.baseline)
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
    args = parser.parse_args()
    if args.blas_interfaces:
        return blas_interfaces(args)
    if args.apply:
        if not args.module:
            raise SystemExit("--apply needs at least one --module")
        return apply_staged(args)
    if args.extract:
        return extract(args)
    return build(args)


if __name__ == "__main__":
    sys.exit(main() or 0)

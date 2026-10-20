"""Track 0: verify that the ctypes layouts of the Python wrapper match the Fortran derived types exactly.

The script imports camb from a built tree (default: the frozen A0 build), walks the ctypes _fields_ of
    EFTCAMB_timestep_cache, EFTCAMB_parameter_cache, EFTCAMB (TEFTCAMB), CAMBparams (with nested
    TransferParams, AccuracyParams, SourceTermParams, CustomSources) and the CAMBdata prefix mapped in python,
generates a Fortran program that prints, for every mapped component, the byte offset (loc(x%comp) - loc(x)),
the storage size and the Fortran type/kind (via select type on class(*)), compiles it against the tree's module
files, links it against the tree's camblib.so, runs it and compares with ctypes offsetof/sizeof.

Allocatable components (array descriptors, class(...), allocatable) are located by allocating them and scanning
the parent object's memory for the data address (gfortran puts the base address first in a descriptor).

Usage: python layout_check.py [--tree PATH] [--build PATH] [--out DIR]
"""

from __future__ import annotations

import argparse
import ctypes
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)

p = argparse.ArgumentParser()
p.add_argument("--tree", default=str(ROOT / "trees" / "A0"))
p.add_argument("--build", default=str(ROOT / "build" / "A0"))
p.add_argument("--out", default=str(ROOT / "work" / "t0-layout" / "layout"))
args = p.parse_args()

TREE, BUILD, OUT = Path(args.tree), Path(args.build), Path(args.out)
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(TREE))
import camb  # noqa: E402
from camb import baseconfig, eftcamb, model, results  # noqa: E402

assert Path(camb.__file__).resolve().is_relative_to(TREE.resolve()), camb.__file__

# python field name -> fortran component name, where they differ (case-insensitive match otherwise)
RENAME = {
    ("CAMBparams", "EFTCAMB_parameter_cache"): "eft_par_cache",
    ("CAMBdata", "Params"): "CP",
    ("CAMBdata", "OnlyTransfers"): "OnlyTransfer",
}
# python structure class -> fortran type name (for nested non-F2003 structures)
FTYPE = {
    "TransferParams": "TransferParams",
    "AccuracyParams": "AccuracyParams",
    "SourceTermParams": "SourceTermParams",
    "CustomSources": "TCustomSourceParams",
}
# python base classes of AllocatableObject fields whose _fortran_class_name_ defaults to the python name
FCLASS = {
    "InitialPower": "TInitialPower",
    "RecombinationModel": "TRecombinationModel",
    "ReionizationModel": "TReionizationModel",
    "DarkEnergyModel": "TDarkEnergyModel",
    "NonLinearModel": "TNonLinearModel",
}


def ftype_of(cls) -> str:
    if cls.__name__ in FCLASS:
        return FCLASS[cls.__name__]
    return getattr(cls, "_fortran_class_name_", None) or FTYPE[cls.__name__]


def leaves(cls, fpath: str, pyoff: int, pypath: str):
    """Yield one record per mapped ctypes leaf field of cls (recursing into nested structures)."""
    for name, ctype in cls._fields_:
        pname = name[1:] if name.startswith("_") else name
        desc = cls.__dict__.get(pname)
        field = getattr(cls, name)
        fname = RENAME.get((cls.__name__, pname), pname)
        rec = {
            "py": f"{pypath}.{pname}",
            "f": f"{fpath}%{fname}",
            "py_off": pyoff + field.offset,
            "py_size": field.size,
        }
        if isinstance(desc, baseconfig.BoolField):
            rec["kind"] = "logical4"
            rec["py_type"] = "c_bool->c_int"
        elif isinstance(ctype, type) and issubclass(ctype, baseconfig.CAMB_Structure):
            # nested structure: record itself, then recurse
            rec["kind"] = f"type:{ftype_of(ctype)}"
            rec["py_type"] = ctype.__name__
            rec["struct"] = True
            yield rec
            yield from leaves(ctype, rec["f"], rec["py_off"], rec["py"])
            continue
        elif isinstance(ctype, type) and issubclass(ctype, baseconfig._AllocatableObjectArray):
            rec["kind"] = "allocobjarray"
            rec["py_type"] = f"AllocatableObjectArray({ctype._baseclass.__name__})"
        elif isinstance(ctype, type) and issubclass(ctype, baseconfig._AllocatableObject):
            rec["kind"] = f"allocobj:{ftype_of(ctype._baseclass)}"
            rec["py_type"] = f"AllocatableObject({ctype._baseclass.__name__})"
        elif ctype is baseconfig.AllocatableArrayDouble:
            rec["kind"] = "allocarr:real8"
            rec["py_type"] = "AllocatableArrayDouble"
        elif ctype is baseconfig.AllocatableArrayInt:
            rec["kind"] = "allocarr:int4"
            rec["py_type"] = "AllocatableArrayInt"
        elif isinstance(ctype, type) and issubclass(ctype, ctypes.Array):
            base = {ctypes.c_double: "real8", ctypes.c_int: "int4"}[ctype._type_]
            rec["kind"] = f"array:{base}:{ctype._length_}"
            rec["py_type"] = f"{ctype._type_.__name__}*{ctype._length_}"
        elif ctype is ctypes.c_double:
            rec["kind"], rec["py_type"] = "real8", "c_double"
        elif ctype is ctypes.c_int:
            named = isinstance(desc, baseconfig.NamedIntField)
            rec["kind"], rec["py_type"] = "int4", "c_int(named)" if named else "c_int"
        elif ctype is ctypes.c_void_p:
            rec["kind"], rec["py_type"] = "funptr", "c_void_p"
        else:
            raise TypeError(f"unhandled ctype {ctype} for {pname}")
        yield rec


TOP = [
    # (label, python class, fortran type, fortran variable)
    ("timestep_cache", eftcamb.EFTCAMB_timestep_cache, "TEFTCAMB_timestep_cache", "lc_tc"),
    ("parameter_cache", eftcamb.EFTCAMB_parameter_cache, "TEFTCAMB_parameter_cache", "lc_pc"),
    ("TEFTCAMB", eftcamb.EFTCAMB, "TEFTCAMB", "lc_ef"),
    ("CAMBparams", model.CAMBparams, "CAMBparams", "lc_cp"),
    ("CAMBdata", results.CAMBdata, "CAMBdata", "lc_cd"),
    ("TransferParams", model.TransferParams, "TransferParams", "lc_tp"),
    ("AccuracyParams", model.AccuracyParams, "AccuracyParams", "lc_ap"),
    ("SourceTermParams", model.SourceTermParams, "SourceTermParams", "lc_sp"),
    ("CustomSources", model.CustomSources, "TCustomSourceParams", "lc_cs"),
]
# Fortran components that python does not map, printed for information (offsets of the unmapped tail)
EXTRA = {
    "TEFTCAMB": [("outroot", "allocchar"), ("model", "allocobj:EFTCAMB_model"),
                 ("read_param_keys", "allocchararr"), ("read_param_values", "allocchararr")],
}

records: dict[str, list[dict]] = {}
for label, cls, ftype, var in TOP:
    records[label] = list(leaves(cls, var, 0, label))

# ---------------------------------------------------------------------------------------------
# generate the Fortran program

KIND_TYPES = [
    ("integer(4)", "int4"), ("integer(8)", "int8"), ("real(4)", "real4"), ("real(8)", "real8"),
    ("logical(4)", "logical4"), ("logical(1)", "logical1"), ("logical(8)", "logical8"),
]
DERIVED = ["TEFTCAMB", "TEFTCAMB_parameter_cache", "TEFTCAMB_timestep_cache", "CAMBparams", "CAMBdata",
           "TransferParams", "AccuracyParams", "SourceTermParams", "TCustomSourceParams"]
CLASSES = ["TInitialPower", "TRecombinationModel", "TReionizationModel", "TDarkEnergyModel", "TNonLinearModel",
           "TSourceWindow", "EFTCAMB_model"]

lines = [
    "program layout_check",
    "use iso_c_binding",
    "use precision",
    "use classes",
    "use SourceWindows",
    "use DarkEnergyInterface",
    "use model",
    "use results",
    "use EFTCAMB_cache",
    "use EFTCAMB_abstract_model",
    "use EFTCAMB_main",
    "use EFTCAMB_pure_EFT_std",
    "implicit none",
    "type(TEFTCAMB_timestep_cache), target, save :: lc_tc",
    "type(TEFTCAMB_parameter_cache), target, save :: lc_pc",
    "type(TEFTCAMB), target, save :: lc_ef",
    "type(CAMBparams), target, save :: lc_cp",
    "type(CAMBdata), pointer :: lc_cd",
    "type(TransferParams), target, save :: lc_tp",
    "type(AccuracyParams), target, save :: lc_ap",
    "type(SourceTermParams), target, save :: lc_sp",
    "type(TCustomSourceParams), target, save :: lc_cs",
    "type(TEFTCAMB_timestep_cache), target, save :: lc_tc2(3)",
    "allocate(lc_cd)",
]


def size_line(label: str, var: str, ftype: str) -> str:
    return f"write(*,'(a,1x,a,1x,i0)') 'SIZE', '{label}', storage_size({var})/8"


def alloc_stmt(rec: dict) -> list[str]:
    k, f = rec["kind"], rec["f"]
    if k.startswith("allocarr"):
        return [f"if (allocated({f})) deallocate({f})", f"allocate({f}(3))"]
    if k == "allocobjarray":
        return [f"if (allocated({f})) deallocate({f})", f"allocate({f}(2))"]
    if k.startswith("allocobj"):
        t = k.split(":", 1)[1]
        if t == "EFTCAMB_model":
            return [f"if (.not. allocated({f})) allocate(EFTCAMB_std_pure_EFT :: {f})"]
        return [f"if (allocated({f})) deallocate({f})", f"allocate({f})"]
    if k == "allocchar":
        return [f"{f} = 'abcdefgh'"]
    if k == "allocchararr":
        return [f"if (allocated({f})) deallocate({f})", f"allocate({f}(3))"]
    return []


def emit(label: str, var: str, rec: dict) -> list[str]:
    k, f = rec["kind"], rec["f"]
    base = f"loc({var})"
    key = rec.get("py", f)
    out = []
    if k.startswith("alloc"):
        out += alloc_stmt(rec)
        if k == "allocchar":
            data = f"loc({f})"
        elif k in ("allocobjarray", "allocchararr") or k.startswith("allocarr"):
            data = f"loc({f}(1))"
        else:
            data = f"loc({f})"
        out.append(f"call scan('{label}', '{key}', {base}, storage_size({var})/8, {data})")
        if k.startswith("allocarr"):
            out.append(f"call kindof('{label}', '{key}#elem', {f}(1))")
        elif k.startswith("allocobj:"):
            out.append(f"call kindof('{label}', '{key}#dyn', {f})")
        return out
    if k.startswith("array"):
        n = int(k.split(":")[2])
        out.append(
            f"write(*,'(a,1x,a,1x,a,1x,i0,1x,i0,1x,i0)') 'OFF', '{label}', '{key}', {f'loc({f}(1))'} - {base}, "
            f"storage_size({f})/8*size({f}), size({f})"
        )
        out.append(f"call kindof('{label}', '{key}', {f}(1))")
        return out
    out.append(
        f"write(*,'(a,1x,a,1x,a,1x,i0,1x,i0,1x,i0)') 'OFF', '{label}', '{key}', loc({f}) - {base}, "
        f"storage_size({f})/8, 1"
    )
    if k == "funptr":
        # c_funptr cannot appear in select type: a specific dummy of type(c_funptr) checks the type at compile time
        out.append(f"call is_funptr('{label}', '{key}', {f})")
    else:
        out.append(f"call kindof('{label}', '{key}', {f})")
    return out


for label, cls, ftype, var in TOP:
    vref = var
    lines.append(size_line(label, vref, ftype))
    for rec in records[label]:
        lines += emit(label, vref, rec)
    for comp, kind in EXTRA.get(label, []):
        rec = {"f": f"{vref}%{comp}", "kind": kind, "py": f"{label}.{comp}[unmapped]"}
        lines += emit(label, vref, rec)
# array stride of the timestep cache as passed to fortran explicit-shape dummies
lines.append("write(*,'(a,1x,i0)') 'STRIDE timestep_cache', loc(lc_tc2(2)) - loc(lc_tc2(1))")
lines.append("contains")
# scan subroutine
lines += [
    "subroutine scan(label, key, base, nbytes, target_addr)",
    "character(*), intent(in) :: label, key",
    "integer(c_intptr_t), intent(in) :: base, target_addr",
    "integer, intent(in) :: nbytes",
    "integer(c_intptr_t), pointer :: words(:)",
    "integer :: i, nfound, first",
    "call c_f_pointer(transfer(base, c_null_ptr), words, [nbytes/8])",
    "nfound = 0",
    "first = -1",
    "do i = 1, nbytes/8",
    "  if (words(i) == target_addr) then",
    "    nfound = nfound + 1",
    "    if (first < 0) first = (i-1)*8",
    "  end if",
    "end do",
    "write(*,'(a,1x,a,1x,a,1x,i0,1x,i0)') 'ALLOC', label, key, first, nfound",
    "end subroutine scan",
    "subroutine is_funptr(label, key, x)",
    "character(*), intent(in) :: label, key",
    "type(c_funptr), intent(in) :: x",
    "write(*,'(a,1x,a,1x,a,1x,a)') 'KIND', label, key, 'funptr'",
    "end subroutine is_funptr",
    "subroutine kindof(label, key, x)",
    "character(*), intent(in) :: label, key",
    "class(*), intent(in) :: x",
    "character(len=64) :: code",
    "select type (x)",
]
for decl, code in KIND_TYPES:
    lines += [f"type is ({decl})", f"  code = '{code}'"]
lines += ["type is (character(*))", "  code = 'char'"]
for t in DERIVED:
    lines += [f"type is ({t})", f"  code = 'type:{t}'"]
for t in CLASSES:
    lines += [f"class is ({t})", f"  code = 'class:{t}'"]
lines += [
    "class default",
    "  code = 'unknown'",
    "end select",
    "write(*,'(a,1x,a,1x,a,1x,a)') 'KIND', label, key, trim(code)",
    "end subroutine kindof",
    "end program layout_check",
]

src = OUT / "layout_check.f90"
src.write_text("\n".join(lines) + "\n")

mods = BUILD / "camb" / "Releaselib"
fmods = BUILD / "forutils" / "Releaselib"
exe = OUT / "layout_check"
cmd = [
    "gfortran", "-cpp", "-ffree-line-length-none", "-fopenmp", "-O0", "-g",
    f"-I{mods}", f"-I{fmods}", str(src), "-o", str(exe),
    str(TREE / "camb" / "camblib.so"), "-Wl,-rpath,/Users/marco/Software/miniforge3/lib",
]
comp = subprocess.run(cmd, capture_output=True, text=True)
(OUT / "compile.log").write_text(" ".join(cmd) + "\n" + comp.stdout + comp.stderr)
if comp.returncode:
    print(comp.stdout + comp.stderr)
    sys.exit("compile failed")
run = subprocess.run([str(exe)], capture_output=True, text=True, cwd=OUT)
(OUT / "fortran_output.txt").write_text(run.stdout + run.stderr)
if run.returncode:
    print(run.stdout[-3000:], run.stderr[-3000:])
    sys.exit("run failed")

# ---------------------------------------------------------------------------------------------
# parse and compare
f_off, f_size, f_kind, f_alloc, f_struct = {}, {}, {}, {}, {}
stride = None
for line in run.stdout.splitlines():
    tok = line.split()
    if not tok:
        continue
    if tok[0] == "SIZE":
        f_struct[tok[1]] = int(tok[2])
    elif tok[0] == "OFF":
        f_off[(tok[1], tok[2])] = (int(tok[3]), int(tok[4]), int(tok[5]))
    elif tok[0] == "KIND":
        f_kind[(tok[1], tok[2])] = tok[3]
    elif tok[0] == "ALLOC":
        f_alloc[(tok[1], tok[2])] = (int(tok[3]), int(tok[4]))
    elif tok[0] == "STRIDE":
        stride = int(tok[2])

alloc_sizes = {
    "allocobj": baseconfig._f_allocatable_size,
    "allocarr": baseconfig._f_allocatable_array_size,
    "allocobjarray": baseconfig._f_allocatable_object_array_size,
}
mismatches = []
report = {"structs": {}, "fields": {}, "stride_timestep_cache": stride,
          "python_alloc_sizes": alloc_sizes}
for label, cls, ftype, var in TOP:
    py_size = ctypes.sizeof(cls)
    fs = f_struct[label]
    exact = label in ("timestep_cache", "parameter_cache", "CAMBparams", "TransferParams", "AccuracyParams",
                      "SourceTermParams", "CustomSources")
    ok = (py_size == fs) if exact else (py_size <= fs)
    report["structs"][label] = {"fortran_type": ftype, "fortran_bytes": fs, "ctypes_bytes": py_size,
                                "rule": "equal" if exact else "prefix (<=)", "ok": ok}
    if not ok:
        mismatches.append(f"{label}: sizeof ctypes {py_size} vs storage_size {fs}")
    rows = []
    for rec in records[label]:
        key = rec["py"]
        k = rec["kind"]
        row = {"py": key, "fortran": rec["f"], "py_type": rec["py_type"], "kind": k, "py_off": rec["py_off"],
               "py_size": rec["py_size"]}
        if k.startswith("alloc"):
            off, nfound = f_alloc[(label, key)]
            row["f_off"] = off
            row["nfound"] = nfound
            fam = "allocobjarray" if k == "allocobjarray" else k.split(":")[0]
            row["ok"] = off == rec["py_off"] and nfound == 1 and rec["py_size"] == alloc_sizes[fam]
            if k.startswith("allocarr"):
                row["f_kind"] = f_kind.get((label, key + "#elem"))
                row["ok"] &= row["f_kind"] == k.split(":")[1]
            else:
                row["f_kind"] = f_kind.get((label, key + "#dyn"))
                if k.startswith("allocobj:"):
                    row["ok"] &= row["f_kind"] in (f"type:{k.split(':')[1]}", f"class:{k.split(':')[1]}")
        else:
            off, size, n = f_off[(label, key)]
            row["f_off"], row["f_size"] = off, size
            row["f_kind"] = f_kind.get((label, key))
            if k.startswith("array"):
                _, base, length = k.split(":")
                row["ok"] = off == rec["py_off"] and size == rec["py_size"] and n == int(length) and row[
                    "f_kind"] == base
            elif k.startswith("type:"):
                row["ok"] = off == rec["py_off"] and size == rec["py_size"] and row["f_kind"] == k
            else:
                row["ok"] = off == rec["py_off"] and size == rec["py_size"] and row["f_kind"] == k
        if not row["ok"]:
            mismatches.append(f"{label}: {row}")
        rows.append(row)
    for comp, kind in EXTRA.get(label, []):
        key = f"{label}.{comp}[unmapped]"
        off, nfound = f_alloc[(label, key)]
        rows.append({"py": key, "fortran": f"{var}%{comp}", "kind": kind, "f_off": off, "nfound": nfound,
                     "info_only": True})
    report["fields"][label] = rows
if stride != ctypes.sizeof(eftcamb.EFTCAMB_timestep_cache):
    mismatches.append(f"timestep_cache array stride {stride} vs ctypes {ctypes.sizeof(eftcamb.EFTCAMB_timestep_cache)}")
report["mismatches"] = mismatches
report["n_fields_checked"] = sum(len(v) for v in records.values())
(OUT / "layout_report.json").write_text(json.dumps(report, indent=1))
for label, s in report["structs"].items():
    print(f"{label:18s} fortran={s['fortran_bytes']:7d} ctypes={s['ctypes_bytes']:7d} rule={s['rule']:12s} ok={s['ok']}")
print("timestep_cache array stride:", stride)
print("python allocatable sizes:", alloc_sizes)
print("fields checked:", report["n_fields_checked"], "mismatches:", len(mismatches))
for m in mismatches:
    print("MISMATCH", m)
for comp in report["fields"]["TEFTCAMB"]:
    if comp.get("info_only"):
        print("unmapped TEFTCAMB tail:", comp)

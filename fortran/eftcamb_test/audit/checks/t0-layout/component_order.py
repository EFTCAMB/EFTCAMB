"""Track 0: compare the ordered list of Fortran components (parsed from the type definitions in the source) with the
python ctypes _fields_ names, so that no Fortran component hidden in an alignment gap goes unnoticed.

For types python allocates or embeds (timestep cache, parameter cache, CAMBparams, nested parameter types) the lists
must be identical; for TEFTCAMB and CAMBdata, python must map an exact prefix of the Fortran components.

Usage: python component_order.py [--tree PATH]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
p = argparse.ArgumentParser()
p.add_argument("--tree", default=str(ROOT / "trees" / "A0"))
p.add_argument("--out", default=str(ROOT / "work" / "t0-layout" / "layout" / "component_order.json"))
args = p.parse_args()
TREE = Path(args.tree)
sys.path.insert(0, str(TREE))
from camb import eftcamb, model, results  # noqa: E402

F = TREE / "fortran"
TYPES = [
    # label, source file, fortran type, python class, rule
    ("timestep_cache", F / "eftcamb" / "03_abstract_EFTCAMB_cache.f90", "TEFTCAMB_timestep_cache",
     eftcamb.EFTCAMB_timestep_cache, "equal"),
    ("parameter_cache", F / "eftcamb" / "03_abstract_EFTCAMB_cache.f90", "TEFTCAMB_parameter_cache",
     eftcamb.EFTCAMB_parameter_cache, "equal"),
    ("TEFTCAMB", F / "eftcamb" / "09_EFTCAMB_main.f90", "TEFTCAMB", eftcamb.EFTCAMB, "prefix"),
    ("CAMBparams", F / "model.f90", "CAMBParams", model.CAMBparams, "equal"),
    ("TransferParams", F / "model.f90", "TransferParams", model.TransferParams, "equal"),
    ("AccuracyParams", F / "model.f90", "AccuracyParams", model.AccuracyParams, "equal"),
    ("CustomSources", F / "model.f90", "TCustomSourceParams", model.CustomSources, "equal"),
    ("SourceTermParams", F / "SourceWindows.f90", "SourceTermParams", model.SourceTermParams, "equal"),
    ("CAMBdata", F / "results.f90", "CAMBdata", results.CAMBdata, "prefix"),
]
RENAME = {"eft_par_cache": "EFTCAMB_parameter_cache", "CP": "Params", "OnlyTransfer": "OnlyTransfers"}


def fortran_components(path: Path, tname: str) -> list[str]:
    text = path.read_text()
    # join continuation lines
    text = re.sub(r"&\s*\n\s*&?", " ", text)
    lines = text.splitlines()
    start = None
    pat = re.compile(rf"^\s*type\b.*?(::)?\s*{tname}\s*$", re.I)
    for i, line in enumerate(lines):
        code = line.split("!", 1)[0]
        if pat.match(code):
            start = i
            break
    if start is None:
        raise ValueError(f"type {tname} not found in {path}")
    comps = []
    for line in lines[start + 1 :]:
        code = line.split("!", 1)[0].strip()
        if not code:
            continue
        low = code.lower()
        if low.startswith("contains") or low.startswith("end type"):
            break
        if low.startswith("#"):
            continue
        if "::" in code:
            decl = code.split("::", 1)[1]
        else:
            m = re.match(r"^(real\s*\(\s*dl\s*\)|real\(dl\)|integer|logical)\s+(.*)$", code, re.I)
            if not m:
                raise ValueError(f"cannot parse component line in {tname}: {code}")
            decl = m.group(2)
        # split on commas not inside parentheses
        depth, cur, parts = 0, "", []
        for ch in decl:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            if ch == "," and depth == 0:
                parts.append(cur)
                cur = ""
            else:
                cur += ch
        parts.append(cur)
        for part in parts:
            name = re.split(r"[=(\s]", part.strip(), maxsplit=1)[0]
            if name:
                comps.append(name)
    return comps


report, failures = {}, []
for label, src, tname, cls, rule in TYPES:
    fcomps = [RENAME.get(c, c) for c in fortran_components(src, tname)]
    pyf = [n[1:] if n.startswith("_") else n for n, _ in cls._fields_]
    fl, pl = [c.lower() for c in fcomps], [c.lower() for c in pyf]
    if rule == "equal":
        ok = fl == pl
    else:
        ok = fl[: len(pl)] == pl
    first_diff = next((i for i, (a, b) in enumerate(zip(fl, pl)) if a != b), None)
    report[label] = {"rule": rule, "n_fortran": len(fl), "n_python": len(pl), "ok": ok, "first_diff": first_diff,
                     "fortran_unmapped_tail": fcomps[len(pl):] if rule == "prefix" else []}
    if first_diff is not None:
        report[label]["diff_context"] = {"fortran": fcomps[first_diff: first_diff + 3],
                                         "python": pyf[first_diff: first_diff + 3]}
    if not ok:
        failures.append(label)
    print(f"{label:18s} fortran={len(fl):4d} python={len(pl):4d} rule={rule:6s} ok={ok}"
          + ("" if first_diff is None else f" first_diff@{first_diff}: {report[label]['diff_context']}"))
    if rule == "prefix":
        tail = report[label]["fortran_unmapped_tail"]
        print(f"    unmapped fortran tail ({len(tail)}): {tail[:12]}{' ...' if len(tail) > 12 else ''}")
Path(args.out).write_text(json.dumps(report, indent=1))
print("FAILURES:", failures)

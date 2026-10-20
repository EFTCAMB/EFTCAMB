"""Check that the EFTCAMB patch applied on CAMB v2 (camb..eftcamb-audit) is equivalent to the pre-rebase patch
(1db4141..eftcamb-pre-rebase), modulo formatting.

Fortran sources are reduced to normalised logical statements (continuations joined, comments dropped, lower case,
no whitespace, '::' removed, old-style relational operators mapped).  For a CAMB file f:

    added_old   = stmts(76b5698:f) - stmts(1db4141:f)      (multiset difference)
    added_new   = stmts(5bef49a:f) - stmts(1d1d3d6:f)
    new_only    = added_new - added_old   -> EFT statements that appear only after the rebase
    old_only    = added_old - added_new   -> EFT statements of the old patch missing after the rebase
    (and the same for statements the patch removed from upstream)

For EFTCAMB's own files (fortran/eftcamb/*.f90) the statements of B1 and A0 are compared directly.

Usage: python patch_equivalence.py [--files a.f90,b.f90] [--eftcamb] > report.txt
"""

from __future__ import annotations

import argparse
import re
import subprocess
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
OLD_BASE, OLD_EFT, NEW_BASE, NEW_EFT = "1db4141", "eftcamb-pre-rebase", "camb", "eftcamb-audit"
OPS = {".eq.": "==", ".ne.": "/=", ".lt.": "<", ".gt.": ">", ".le.": "<=", ".ge.": ">="}


def show(rev: str, path: str) -> str:
    r = subprocess.run(["git", "-C", str(REPO), "show", f"{rev}:{path}"], capture_output=True, text=True,
                       errors="replace")
    return r.stdout if r.returncode == 0 else ""


def strip_comment(line: str) -> str:
    if line.lstrip().startswith("!$"):
        return line
    out, q = [], None
    for ch in line:
        if q:
            out.append(ch)
            if ch == q:
                q = None
        elif ch in ("'", '"'):
            q = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out)


def statements(src: str) -> list[str]:
    stmts, buf = [], ""
    for raw in src.splitlines():
        if raw.lstrip().startswith("#"):
            stmts.append(re.sub(r"\s+", "", raw.lower()))
            continue
        line = strip_comment(raw).rstrip()
        if not line.strip():
            continue
        s = line.strip()
        if s.startswith("&"):
            s = s[1:]
        if s.endswith("&"):
            buf += s[:-1]
            continue
        buf += s
        for part in buf.split(";"):
            t = re.sub(r"\s+", "", part.lower()).replace("::", "")
            for k, v in OPS.items():
                t = t.replace(k, v)
            if t:
                stmts.append(t)
        buf = ""
    return stmts


def diff(a: list[str], b: list[str]) -> tuple[Counter, Counter]:
    ca, cb = Counter(a), Counter(b)
    return cb - ca, ca - cb  # added (in b), removed (from a)


def report_camb(path: str) -> None:
    ob, oe, nb, ne = (statements(show(r, path)) for r in (OLD_BASE, OLD_EFT, NEW_BASE, NEW_EFT))
    add_old, rem_old = diff(ob, oe)
    add_new, rem_new = diff(nb, ne)
    sections = [
        ("EFT-added only after rebase", add_new - add_old),
        ("EFT-added only before rebase (dropped?)", add_old - add_new),
        ("upstream stmts removed by EFT only after rebase", rem_new - rem_old),
        ("upstream stmts removed by EFT only before rebase", rem_old - rem_new),
    ]
    print(f"==== {path}: added old {sum(add_old.values())} new {sum(add_new.values())}; "
          f"removed old {sum(rem_old.values())} new {sum(rem_new.values())}")
    for title, c in sections:
        if c:
            print(f"  -- {title}: {sum(c.values())}")
            for s, n in sorted(c.items()):
                print(f"     [{n}] {s[:300]}")


def report_eft(path: str) -> None:
    a, b = statements(show(OLD_EFT, path)), statements(show(NEW_EFT, path))
    add, rem = diff(a, b)
    print(f"==== {path}: B1 {len(a)} stmts, A0 {len(b)} stmts, only-A0 {sum(add.values())}, only-B1 {sum(rem.values())}")
    for title, c in (("only in A0", add), ("only in B1", rem)):
        if c:
            print(f"  -- {title}: {sum(c.values())}")
            for s, n in sorted(c.items()):
                print(f"     [{n}] {s[:300]}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--files", default="")
    p.add_argument("--eftcamb", action="store_true")
    a = p.parse_args()
    if a.eftcamb:
        r = subprocess.run(["git", "-C", str(REPO), "ls-tree", "-r", "--name-only", NEW_EFT, "fortran/eftcamb"],
                           capture_output=True, text=True)
        files = [f for f in r.stdout.split() if f.endswith(".f90")]
        if a.files:
            files = [f for f in files if Path(f).name in a.files.split(",")]
        for f in files:
            report_eft(f)
    else:
        files = a.files.split(",") if a.files else [
            "camb.f90", "cmbmain.f90", "halofit.f90", "equations.f90", "massive_neutrinos.f90", "lensing.f90",
            "reionization.f90", "model.f90", "results.f90", "camb_python.f90", "DarkEnergyInterface.f90",
            "SourceWindows.f90", "subroutines.f90", "config.f90", "recfast.f90", "transfer.f90", "Makefile_main",
            "Makefile"]
        for f in files:
            report_camb(f"fortran/{f}")


if __name__ == "__main__":
    main()

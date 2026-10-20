"""Compare the 'EFTCAMB MOD START ... END' blocks in the CAMB sources of two trees (A0 vs B1).

Blocks are normalised (lower case, all whitespace removed, '&' continuations joined, comments dropped) so that the
forformat reformatting does not show up.  Blocks are matched by their normalised content; blocks present in only one
tree are printed with file and line number.

Usage: python mod_blocks_diff.py TREE_A TREE_B
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

FILES = ["camb.f90", "cmbmain.f90", "halofit.f90", "equations.f90", "massive_neutrinos.f90", "lensing.f90",
         "reionization.f90", "model.f90", "results.f90", "camb_python.f90", "DarkEnergyInterface.f90",
         "SourceWindows.f90", "subroutines.f90", "config.f90", "recfast.f90", "transfer.f90", "bessels.f90"]


def norm_line(s: str) -> str:
    s = s.split("!", 1)[0] if not s.lstrip().startswith("!$") else s
    s = s.lower().replace("&", "")
    return re.sub(r"\s+", "", s)


def blocks(path: Path) -> list[tuple[int, str, str]]:
    out = []
    if not path.exists():
        return out
    lines = path.read_text(errors="replace").splitlines()
    i = 0
    while i < len(lines):
        if "EFTCAMB MOD START" in lines[i]:
            j = i + 1
            while j < len(lines) and "EFTCAMB MOD END" not in lines[j]:
                j += 1
            raw = "\n".join(lines[i:j + 1])
            normd = "".join(norm_line(x) for x in lines[i:j + 1])
            out.append((i + 1, normd, raw))
            i = j + 1
        else:
            i += 1
    return out


def main() -> None:
    ta, tb = Path(sys.argv[1]), Path(sys.argv[2])
    for f in FILES:
        ba, bb = blocks(ta / "fortran" / f), blocks(tb / "fortran" / f)
        if not ba and not bb:
            continue
        ca, cb = Counter(x[1] for x in ba), Counter(x[1] for x in bb)
        only_a = [x for x in ba if ca[x[1]] > cb.get(x[1], 0)]
        only_b = [x for x in bb if cb[x[1]] > ca.get(x[1], 0)]
        print(f"==== {f}: A blocks {len(ba)}, B blocks {len(bb)}, only-A {len(only_a)}, only-B {len(only_b)}")
        for ln, _, raw in only_a:
            print(f"--- only in A {f}:{ln}\n{raw}\n")
        for ln, _, raw in only_b:
            print(f"--- only in B {f}:{ln}\n{raw}\n")


if __name__ == "__main__":
    main()

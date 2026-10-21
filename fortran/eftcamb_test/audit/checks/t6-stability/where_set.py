"""For each EFTCAMB model source, report in which subroutine/function each cache field is assigned."""
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
SRC = REPO / "fortran" / "eftcamb"
fields = sys.argv[1:] or ["EFTGamma2PPP", "EFTGamma3PPP", "EFTGamma3PPPP", "EFTLambdadotdot", "Hdotdotdot",
                          "EFTOmegaPPPP", "EFTGamma2PP", "EFTGamma3PP", "EFTOmegaPPP", "EFTcdot", "Hdotdot"]
start = re.compile(r"^\s*(?:recursive\s+|pure\s+|elemental\s+)*(subroutine|function)\s+(\w+)", re.I)
endre = re.compile(r"^\s*end\s+(subroutine|function)", re.I)
for f in sorted(SRC.rglob("*.f90")):
    cur = None
    found = {}
    for line in f.read_text(errors="replace").splitlines():
        code = line.split("!", 1)[0]
        m = start.match(code)
        if m:
            cur = m.group(2)
            continue
        if endre.match(code):
            cur = None
            continue
        for fld in fields:
            if re.search(rf"eft_cache%{fld}\s*=(?!=)", code):
                found.setdefault(fld, set()).add(cur)
    if found:
        print(f.relative_to(SRC))
        for fld in fields:
            if fld in found:
                print(f"    {fld:16s} {sorted(found[fld], key=str)}")

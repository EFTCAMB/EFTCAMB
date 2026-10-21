"""t3-gr (a): list every EFTCAMB MOD block in the CAMB sources (not fortran/eftcamb) with its first code lines, to
classify which ones are active when EFTflag = 0."""
import re
import sys
from pathlib import Path

repo = Path(sys.argv[1])
files = sorted(repo.glob("fortran/*.f90")) + sorted(repo.glob("camb/*.py")) + [repo / "fortran" / "Makefile"]
for f in files:
    lines = f.read_text(errors="replace").splitlines()
    i = 0
    while i < len(lines):
        if "EFTCAMB MOD START" in lines[i]:
            label = lines[i].split("EFTCAMB MOD START", 1)[1].strip(" :.")
            code = []
            j = i + 1
            while j < len(lines) and "EFTCAMB MOD END" not in lines[j]:
                s = lines[j].strip()
                if s and not s.startswith(("!", "#")) or s.startswith(("#ifdef", "#if", "#endif", "#else")):
                    code.append(s)
                j += 1
            print(f"{f.relative_to(repo)}:{i + 1}-{j + 1} [{label[:70]}] :: " + " | ".join(code[:3])[:230])
            i = j
        i += 1

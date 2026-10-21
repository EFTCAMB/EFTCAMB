"""Verify (mass-stability-verdict-unreliable): build t6's fixaddl tree (F-009 candidate fix) plus a dimensionally
consistent EFTcPP = d^2 EFTc/da^2 in EFTCAMBModelComputeMassFactors, to test the EFTcPP root cause in the real code."""

import sys
from pathlib import Path

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
SRC = f"{ROOT}/work/t6-stability/fixaddl"
DEST = Path(f"{ROOT}/work/verify/massverd/fixaddl_cpp")
OLD = ("    EFTcPP = (eft_cache%EFTcdotdot/eft_cache%adotoa + (3._dl - eft_cache%Hdot/eft_cache%adotoa)*eft_cache%EFTcdot + &\n"
       "        2._dl*eft_cache%adotoa*eft_cache%EFTc)/(a**2*eft_cache%adotoa**2)\n")
NEW = ("    EFTcPP = (eft_cache%EFTcdotdot + (3._dl*eft_cache%adotoa - eft_cache%Hdot/eft_cache%adotoa)*eft_cache%EFTcdot + &\n"
       "        2._dl*eft_cache%adotoa**2*eft_cache%EFTc)/(a**2*eft_cache%adotoa**2)\n")
if not DEST.exists():
    import shutil
    shutil.copytree(SRC, DEST, symlinks=True)
f = DEST / "fortran/eftcamb/06_abstract_EFTCAMB_model.f90"
t = f.read_text()
if OLD in t:
    f.write_text(t.replace(OLD, NEW))
    print("patched EFTcPP")
elif NEW in t:
    print("already patched")
else:
    sys.exit("EFTcPP line not found")
print(harness.build_tree(DEST, DEST))

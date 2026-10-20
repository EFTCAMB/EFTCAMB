"""Adversarial verifier (t3-extgal-gr-point-background): build A0 copies with candidate ExtGal fixes.

hddfix: only the Hdotdot denominator ((S+2) H^2 - a^2 h0^2 S Omega_tot)^3 in ComputeHubbleDer.
allfix: hddfix + t2d-static FIX3 (massive nu inside ComputeAdotoa, real bisection tolerance) + tracker at all a
        (the a < 0.2 cache branch disabled).
Patched sources are prepared in work/verify/extgal_gr/src/{hddfix,allfix}.f90.
"""

import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work/verify/extgal_gr"
REL = Path("fortran/eftcamb/08f_full_models/008p7_Extended_Galileon.f90")
for name in sys.argv[1:]:
    dest = W / name
    if not dest.exists():
        shutil.copytree(harness.TREES["A0"], dest, symlinks=True)
    shutil.copy(W / "src" / f"{name}.f90", dest / REL)
    print(name, harness.build_tree("A0", dest))

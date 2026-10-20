"""Verifier (ic1-pidot-transient-noisy-coefficients): recreate work/verify/ic1noisy/A0_icdump = A0 + print-only
instrumentation of 09_EFTCAMB_IC.f90 case(1), by applying work/verify/ic1noisy/A0_icdump_instrumentation.diff to a fresh
copy of the frozen A0 tree, then build it (release flags) with harness.build_tree."""

import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy"
dest = W / "A0_icdump"
if not dest.exists():
    shutil.copytree(harness.tree_path("A0"), dest, symlinks=True)
    (dest / "camb" / "camblib.so").unlink()
    subprocess.run(["patch", str(dest / "fortran" / "eftcamb" / "09_EFTCAMB_IC.f90"),
                    str(W / "A0_icdump_instrumentation.diff")], check=True)
for _ in range(3):  # forutils' parallel make can race on a clean build dir; retry
    r = harness.build_tree("A0", dest)
    if (dest / "camb" / "camblib.so").exists():
        break
print(r)

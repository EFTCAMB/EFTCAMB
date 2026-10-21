"""Make a copy of tree A0 with a minimal root-cause patch for the cubic spline (flag 9) finding and build it.

Patch (verifier's own, not the reporter's FIX.diff):
  1. 04p016_spline_parametrizations_1D.f90: strict end tests `x <= x_initial .or. x >= x_final` -> `x < ... .or. x > ...`
  2. 02_interpolation_linear_1D.f90: drop the 8 lines forcing yp/ypp/yppp/ypppp to 0 at the first/last node
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

DEST = harness.AUDIT_ROOT / "work" / "verify" / "spline_cubic_fixmin"
if not DEST.exists():
    shutil.copytree(harness.TREES["A0"], DEST, symlinks=True)
    f1 = DEST / "fortran/eftcamb/04f_parametrizations_1D/04p016_spline_parametrizations_1D.f90"
    s = f1.read_text()
    old = "x <= self%interpolation%x_initial .or. x >= self%interpolation%x_final"
    assert s.count(old) == 6, s.count(old)
    f1.write_text(s.replace(old, "x < self%interpolation%x_initial .or. x > self%interpolation%x_final"))
    f2 = DEST / "fortran/eftcamb/02_interpolation_linear_1D.f90"
    lines = f2.read_text().splitlines(keepends=True)
    keep = [ln for ln in lines if not any(ln.strip() == f"self%{v}({i}) = 0._dl"
                                          for v in ("yp", "ypp", "yppp", "ypppp") for i in ("1", "self%num_points"))]
    assert len(lines) - len(keep) == 8, len(lines) - len(keep)
    f2.write_text("".join(keep))
print(harness.build_tree("A0", DEST))

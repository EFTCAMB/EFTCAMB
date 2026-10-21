"""Verifier (spline-integral-convention): copy tree A0, replace ONLY the bodies of Spline1DIntegral (04p016) and
SplineQ1DIntegral (04p019) with the verifier's own quadrature (5-point Gauss-Legendre on uniform 0.1 pieces in
N = ln a, no node splitting; independent of the reporter's FIX.diff), and build it.

I(a) = a^2 exp(-3 int_0^{ln a} (1 + w(N)) dN)  (EFTCAMB convention of every analytic w(a) parametrization)

Nothing else is changed (value()/derivatives untouched), so if designer f(R) / shift-symmetric runs with a spline w
then reproduce EFTwDE = 1 the defect is fully attributed to integral().
usage: /Users/marco/Software/miniforge3/bin/python patch_spline_integral_only.py
"""

import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

DEST = harness.AUDIT_ROOT / "work" / "verify" / "spline_integral" / "tree_intfix"
PAR = "fortran/eftcamb/04f_parametrizations_1D/"

BODY = """
    real(dl), parameter :: gx(5) = [0._dl, -0.5384693101056831_dl, 0.5384693101056831_dl, &
        -0.9061798459386640_dl, 0.9061798459386640_dl]
    real(dl), parameter :: gw(5) = [0.5688888888888889_dl, 0.4786286704993665_dl, 0.4786286704993665_dl, &
        0.2369268850561891_dl, 0.2369268850561891_dl]
    real(dl) :: n0, lo, hi, h, s, c
    integer  :: i, j, npc

    {name} = 0._dl
    if (x <= 0._dl) return
    n0 = log(x)
    lo = min(n0, 0._dl)
    hi = max(n0, 0._dl)
    npc = max(1, ceiling((hi - lo)/0.1_dl))
    h = (hi - lo)/npc
    s = 0._dl
    do i = 1, npc
        c = lo + (i - 0.5_dl)*h
        do j = 1, 5
            s = s + gw(j)*0.5_dl*h*(1._dl + self%value(exp(c + 0.5_dl*h*gx(j))))
        end do
    end do
    if (n0 < 0._dl) then
        {name} = exp(2._dl*n0 + 3._dl*s)
    else
        {name} = exp(2._dl*n0 - 3._dl*s)
    end if

    end function {name}
"""


def patch(path: Path, name: str):
    s = path.read_text()
    # keep the declaration block up to and including "real(dl) :: <name> ..." line, replace the rest of the function
    pat = re.compile(r"(    function " + name + r"\(self, x, eft_cache\).*?real\(dl\) :: " + name + r"[^\n]*\n)"
                     r"(.*?)\n    end function " + name + r"\n", re.S)
    m = pat.search(s)
    assert m, name
    s = s[: m.start()] + m.group(1) + BODY.format(name=name) + s[m.end():]
    path.write_text(s)


def main():
    if not DEST.exists():
        shutil.copytree(harness.TREES["A0"], DEST, symlinks=True)
        patch(DEST / PAR / "04p016_spline_parametrizations_1D.f90", "Spline1DIntegral")
        patch(DEST / PAR / "04p019_spline5_parametrizations_1D.f90", "SplineQ1DIntegral")
    print(harness.build_tree("A0", DEST))


if __name__ == "__main__":
    main()

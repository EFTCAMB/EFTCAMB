"""Track t4b-derivs: apply the candidate fix for the Scaling Cubic Galileon third-derivative formulas to a TREE COPY.

Fixes (verified symbolically by scg_formula_sympy.py):
  * hpprime: overall sign (+1/D instead of -1/D) and y2 coefficient (1 - 2 A lambda + 2 A beta2), all three copies;
  * fpprime: add the missing -(2 A lambda - 1)(sqrt6 - 6 A x1) x1pp term inside the -3/D bracket, both copies;
  * EFTLambdadotdot = -2 adotoa^4 ((1 + 3 h) K + K'), with EFTLambdadot = -2 adotoa^3 K;
  * EFTGamma2PPP: 3*h*2 -> 3*h**2.

    python patch_scg_derivs.py TREE_COPY      (never run on the repository checkout)
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]


def patch(text: str) -> str:
    def rep(old: str, new: str, count: int) -> None:
        nonlocal text
        n = text.count(old)
        if n != count:
            raise SystemExit(f"expected {count} occurrences, found {n}: {old[:80]!r}")
        text = text.replace(old, new)

    rep("    hpprime = -1._dl/(1._dl + 2._dl*self%SCG_A*(3._dl*self%SCG_A - self%SCG_lambda))* &",
        "    hpprime = +1._dl/(1._dl + 2._dl*self%SCG_A*(3._dl*self%SCG_A - self%SCG_lambda))* &", 3)
    rep("(1._dl + 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*self%SCG_A*self%SCG_beta2)*(y2p**2._dl + y2*y2pp)))",
        "(1._dl - 2._dl*self%SCG_A*self%SCG_lambda + 2._dl*self%SCG_A*self%SCG_beta2)*(y2p**2._dl + y2*y2pp)))", 3)
    rep("        (-self%SCG_A*OmegaRadpprime + 6._dl*self%SCG_A*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl)*x1p**2._dl + &",
        "        (-self%SCG_A*OmegaRadpprime - &\n"
        "        (2._dl*self%SCG_A*self%SCG_lambda - 1._dl)*(sqrt(6._dl) - 6._dl*self%SCG_A*x1)*x1pp + &\n"
        "        6._dl*self%SCG_A*(2._dl*self%SCG_A*self%SCG_lambda - 1._dl)*x1p**2._dl + &", 2)
    rep("x1p*(2._dl - 6._dl*h + 3._dl*h*2._dl + 3._dl*hprime)", "x1p*(2._dl - 6._dl*h + 3._dl*h**2._dl + 3._dl*hprime)", 1)
    # EFTLambdadotdot: replace the whole statement (from its first line to the end of the subroutine body)
    start = text.index("    eft_cache%EFTLambdadotdot = -2._dl*adotoa**4._dl* &")
    end = text.index("    end subroutine EFTCAMBScalingCubicBackgroundEFTFunctions")
    new_stmt = (
        "    LK = 2._dl*sqrt(6._dl)*self%SCG_A*x1*h**2 + &\n"
        "        3._dl*h*((-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1**2 + y1**2 + y2**2 + sqrt(6._dl)*self%SCG_A*x1p) + &\n"
        "        x1*(sqrt(6._dl)*self%SCG_A*hprime + 3._dl*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1p) + &\n"
        "        3._dl*y1*y1p + 3._dl*y2*y2p + sqrt(6._dl)*self%SCG_A*x1pp\n"
        "    LKp = 2._dl*sqrt(6._dl)*self%SCG_A*(x1p*h**2 + 2._dl*x1*h*hprime) + &\n"
        "        3._dl*hprime*((-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1**2 + y1**2 + y2**2 + &\n"
        "        sqrt(6._dl)*self%SCG_A*x1p) + &\n"
        "        3._dl*h*(2._dl*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1*x1p + 2._dl*y1*y1p + 2._dl*y2*y2p + &\n"
        "        sqrt(6._dl)*self%SCG_A*x1pp) + &\n"
        "        x1p*(sqrt(6._dl)*self%SCG_A*hprime + 3._dl*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1p) + &\n"
        "        x1*(sqrt(6._dl)*self%SCG_A*hpprime + 3._dl*(-1._dl + 2._dl*self%SCG_A*self%SCG_lambda)*x1pp) + &\n"
        "        3._dl*(y1p**2 + y1*y1pp) + 3._dl*(y2p**2 + y2*y2pp) + sqrt(6._dl)*self%SCG_A*x1ppp\n"
        "    ! EFTLambdadot = -2 adotoa^3 LK, EFTLambdadotdot = adotoa (d EFTLambdadot/dN - 2 EFTLambdadot)\n"
        "    eft_cache%EFTLambdadotdot = -2._dl*adotoa**4*((1._dl + 3._dl*h)*LK + LKp)\n"
    )
    text = text[:start] + new_stmt + text[end:]
    decl = ("    real(dl) :: f, fprime, fpprime, h, hprime, hpprime, grhom, grhorad, H2, EFTcdot, EFTLAMBDAdot, EFTc, "
            "EFTLAMBDA, &\n        adotoa, Hdot\n")
    if decl not in text:
        raise SystemExit("declaration block not found")
    text = text.replace(decl, decl + "    real(dl) :: LK, LKp\n", 1)
    return text


def main() -> None:
    tree = Path(sys.argv[1]).resolve()
    if tree == REPO or REPO in tree.parents:
        raise SystemExit("refusing to patch the repository checkout")
    f = tree / "fortran/eftcamb/08f_full_models/008p6_Scaling_Cubic.f90"
    f.write_text(patch(f.read_text()))
    print("patched", f)


if __name__ == "__main__":
    main()

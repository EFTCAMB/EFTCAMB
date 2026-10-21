"""Accuracy of the Extended Galileon background (008p7_Extended_Galileon.f90:479).

EFTCAMBExtendedGalileonComputeAdotoa solves  Omega_phi0 + Omega_tot*(x/a)^S - (x/a)^(2+S) = 0  for x = adotoa/H0 by
bisection on [0, 1e9] with the loop condition  `sqrt(fsoluction**2) > 10**(-1) .and. counter < 50**1`.
10**(-1) is an integer expression and evaluates to 0, so the tolerance is zero and the loop always stops after 50
halvings: the bracket is then 1e9/2^50 ~ 9e-7 wide, which bounds the accuracy of adotoa for a >= 0.2 (for a < 0.2 a
different, closed formula is used).

This script (tree A0, background only, stability checks off) compares the code's adotoa with the exact root of the
same equation (brentq at 1e-15) on a grid of a, and reports the relative error and the jump at a = 0.2.
Usage: python check_extgal_bisection.py [TREE]
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else harness.TREES["A0"]
INI = harness.REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"


def main() -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    import numpy as np
    from scipy.optimize import brentq

    sys.path.insert(0, str(TREE))
    import camb

    out_dir = harness.AUDIT_ROOT / "work/t2d-static/extgal"
    out_dir.mkdir(parents=True, exist_ok=True)
    w = out_dir / "w.ini"
    w.write_text(f"DEFAULT({INI})\nfeedback_level = 0\nEFTCAMB_feedback_level = 0\nEFTCAMB_write_background = F\n"
                 f"EFT_ghost_stability = F\nEFT_gradient_stability = F\nEFT_additional_priors = F\n"
                 f"output_root = {out_dir / 'run'}\n")
    with harness.slot(1):
        pars = camb.read_ini(str(w), no_validate=True)
        res = camb.get_background(pars)
        a = np.concatenate([np.linspace(0.2, 1.0, 801), [0.2 - 1e-9, 0.2 + 1e-9]])
        _, arr = pars.EFTCAMB.get_eft_functions(res, a)
        pc = res.Params.EFTCAMB_parameter_cache
        S = float(res.Params.EFTCAMB.ExtendedGalileon_S) if hasattr(res.Params.EFTCAMB, "ExtendedGalileon_S") else 2.0
    h0 = pc.h0_Mpc
    om_m = pc.omegac + pc.omegab
    om_r = pc.omegag + pc.omegar
    om_phi = pc.omegav
    grhonu = np.asarray(arr["grhonu_tot"]) if "grhonu_tot" in arr.dtype.names else np.zeros_like(a)
    code = np.asarray(arr["adotoa"])
    exact = np.empty_like(a)
    for i, ai in enumerate(a):
        om_tot = om_m * ai**-3 + om_r * ai**-4 + grhonu[i] / (3 * h0**2 * ai**2)

        def f(x):
            return om_phi + om_tot * (x / ai) ** S - (x / ai) ** (2 + S)

        exact[i] = brentq(f, 1e-6, 1e3, xtol=1e-16, rtol=1e-15) * h0
    rel = np.abs(code - exact) / exact
    main_grid = slice(0, 801)
    print(f"tree {TREE}; S = {S}; omega_phi0 = {om_phi:.6f}; massive-nu term max = {grhonu.max():.3e}")
    print(f"a in [0.2, 1]: max |adotoa_code/adotoa_exact - 1| = {rel[main_grid].max():.3e}, "
          f"median = {np.median(rel[main_grid]):.3e}")
    print(f"bisection bound 1e9/2^51 / (adotoa/H0) at a = 1: {1e9 / 2**51 / (exact[800] / h0):.3e}")
    print(f"jump at a = 0.2 (closed formula below, bisection above): adotoa(0.2-1e-9)/adotoa(0.2+1e-9) - 1 = "
          f"{code[801] / code[802] - 1:+.3e}  (exact: {exact[801] / exact[802] - 1:+.3e})")
    # noise in the numerical derivative of adotoa caused by the quantised bisection result
    dlnH = np.diff(np.log(code[main_grid])) / np.diff(a[main_grid])
    dlnH_ex = np.diff(np.log(exact[main_grid])) / np.diff(a[main_grid])
    print(f"d ln(adotoa)/da from code vs exact root: max rel diff = "
          f"{np.max(np.abs(dlnH - dlnH_ex) / np.abs(dlnH_ex)):.3e}")
    print("expected: relative error at the 1e-10 level or below if the loop converged to a real tolerance")


if __name__ == "__main__":
    main()


def detail() -> None:
    """Print code vs exact adotoa at a few scale factors, with and without the massive-neutrino term."""
    os.environ["OMP_NUM_THREADS"] = "1"
    import numpy as np
    from scipy.optimize import brentq

    sys.path.insert(0, str(TREE))
    import camb

    out_dir = harness.AUDIT_ROOT / "work/t2d-static/extgal"
    w = out_dir / "w.ini"
    with harness.slot(1):
        pars = camb.read_ini(str(w), no_validate=True)
        res = camb.get_background(pars)
        a = np.array([0.1, 0.19, 0.2 - 1e-9, 0.2 + 1e-9, 0.3, 0.5, 0.8, 1.0])
        _, arr = pars.EFTCAMB.get_eft_functions(res, a)
        pc = res.Params.EFTCAMB_parameter_cache
    h0 = pc.h0_Mpc
    print("omegac, omegab, omegag, omegar, omegav, omegan, h0_Mpc:", pc.omegac, pc.omegab, pc.omegag, pc.omegar,
          pc.omegav, pc.omegan, h0)
    for i, ai in enumerate(a):
        for nu in (True, False):
            om_tot = (pc.omegac + pc.omegab) * ai**-3 + (pc.omegag + pc.omegar) * ai**-4 + \
                (arr["grhonu_tot"][i] / (3 * h0**2 * ai**2) if nu else 0.0)
            x = brentq(lambda x: pc.omegav + om_tot * (x / ai) ** 2 - (x / ai) ** 4, 1e-6, 1e3, xtol=1e-16)
            if nu:
                ex_nu = x * h0
            else:
                ex_no = x * h0
        print(f"a={ai:.9f} code={arr['adotoa'][i]:.10e} exact(nu)={ex_nu:.10e} exact(no nu)={ex_no:.10e} "
              f"code/exact-1={arr['adotoa'][i] / ex_nu - 1:+.3e} grhonu_tot={arr['grhonu_tot'][i]:.3e}")


if __name__ == "__main__" and os.environ.get("EXTGAL_DETAIL"):
    detail()

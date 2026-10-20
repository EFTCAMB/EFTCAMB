"""Independent symbolic + numeric check of EFTCAMBFullModelComputeHubbleDer (06p1_abstract_EFTCAMB_full_map.f90:121-165).

1. sympy: build Hdot(tau) exactly as coded, differentiate in conformal time with
       a' = a H,  H' = Hdot,  O^(n)' = a H O^(n+1),
       P' = Pd + 2 H P, Pd' = Pdd + 2 H Pd  (gpresm_t, gpresdotm_t, gpresdotdotm_t: X = x a^2, Xdot = a^2 x'),
       L' = Ld + 2 H L, Ld' = Ldd + 2 H Ld  (EFTLambda, EFTLambdadot, EFTLambdadotdot)
   and compare with the coded Hdotdot and the coded Hdotdotdot (transcribed term by term from the source).
2. numeric: evaluate the coded and the exact Hdotdotdot on the stored cache fields of runs/A0/<model>.npz, with
   EFTLambdadotdot / EFTOmegaPPPP either as stored (0) or replaced by spline derivatives, and compare with the
   spline derivative of the stored Hdotdot.  This splits the error into "missing fields" and "formula" parts.

Usage: python hderivs_hddd_sympy.py
"""

import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.interpolate import CubicSpline

RUNS = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
            "scratchpad/audit/runs/A0")

a, H, Hd, Hdd, O, O1, O2, O3, O4, P, Pd, Pdd, L, Ld, Ldd = sp.symbols(
    "a H Hd Hdd O O1 O2 O3 O4 P Pd Pdd L Ld Ldd", real=True)
O5 = sp.Symbol("O5")


def ddtau(expr, Hdot_val, Hdotdot_val=None):
    """Total conformal-time derivative under the EFTCAMB conventions."""
    rules = {a: a * H, H: Hd, O: a * H * O1, O1: a * H * O2, O2: a * H * O3, O3: a * H * O4, O4: a * H * O5,
             P: Pd + 2 * H * P, Pd: Pdd + 2 * H * Pd, L: Ld + 2 * H * L, Ld: Ldd + 2 * H * Ld}
    res = sum(sp.diff(expr, s) * v for s, v in rules.items())
    res = res + sp.diff(expr, Hd) * (Hdotdot_val if Hdotdot_val is not None else Hdd)
    return res


opo = 1 + O
# coded Hdot (lines 129-131)
Hdot_code = 1 / (opo + a * O1 / 2) * (-(opo + 2 * a * O1 + a**2 * O2) * H**2 / 2 - P / 2 - L / 2)
# coded Hdotdot (lines 133-139)
Hdotdot_code = 1 / (opo + a * O1 / 2) * (
    -a * H**3 * (3 * O1 + 4 * a * O2 + a**2 * O3) / 2
    - H * Hd * (opo + sp.Rational(7, 2) * a * O1 + sp.Rational(3, 2) * a**2 * O2)
    - (Pd + 2 * H * P + Ld + 2 * H * L) / 2)
# coded Hdotdotdot (lines 141-163), transcribed literally
denom = opo + a * O1 / 2
N1 = -a * H * (3 * O1 + a * O2 / 2)
adotoa2 = H**2
Hddd_code = N1 / 2 / denom**2 * (
    -a * H**3 * (3 * O1 + 4 * a * O2 + a**2 * O3) / 2
    - H * Hd * (opo + sp.Rational(7, 2) * a * O1 * sp.Rational(3, 2) * a**2 * O2)
    - (Pdd + 2 * H * Pd - (Ld + 2 * H * L) / 2) / 2) + 1 / denom * (
    -a * adotoa2 * (adotoa2 + 3 * H) * (3 * O1 + 4 * a * O2 + a**2 * O3)
    - a**2 * H**4 * (7 * O2 + 6 * a * O3 + a**2 * O4) / 2
    - (Hd**2 + H * Hdd) * (opo + sp.Rational(7, 2) * a * O1 * sp.Rational(3, 2) * a**2 * O2)
    - a * adotoa2 * Hd * (9 * O1 + 13 * a * O2 + 3 * a**2 * O3) / 2
    - (Pdd + 4 * H * Pd + 2 * Hd * P + 4 * adotoa2 * P) / 2
    - (Ldd + 4 * H * Ld + 2 * Hd * L + 4 * adotoa2 * L) / 2)

# 1a: control  d(Hdot)/dtau == Hdotdot (on shell: Hd -> Hdot_code)
dHdot = ddtau(Hdot_code, Hdot_code)
ctrl = sp.simplify((dHdot - Hdotdot_code).subs(Hd, Hdot_code))
print("sympy control: d(Hdot_code)/dtau - Hdotdot_code (on shell) =", ctrl)
# 1b: exact third derivative = d(Hdotdot_code)/dtau, with Hd, Hdd treated as symbols then put on shell
dHdd = ddtau(Hdotdot_code, Hdot_code, Hdd)
Hddd_exact = dHdd
onshell = {Hdd: Hdotdot_code}
diff_expr = sp.simplify(((Hddd_code - Hddd_exact).subs(onshell)).subs(Hd, Hdot_code))
print("sympy: coded Hdotdotdot - exact d(Hdotdot)/dtau (on shell) == 0 ?", diff_expr == 0)
# at Omega == 0 (quintessence/ADE): is the coded formula exact?
d0 = sp.simplify(diff_expr.subs({O: 0, O1: 0, O2: 0, O3: 0, O4: 0}))
print("sympy: same difference with Omega == 0 identically:", d0)
# does the exact expression depend on O5 (it must not)?
print("sympy: exact expression depends on O5 ?", Hddd_exact.has(O5))

# 1c: which coded pieces are wrong: test the four suspected corrections individually
Hddd_fixed = Hddd_exact
fx = sp.lambdify((a, H, Hd, Hdd, O, O1, O2, O3, O4, P, Pd, Pdd, L, Ld, Ldd), Hddd_code, "numpy")
fe = sp.lambdify((a, H, Hd, Hdd, O, O1, O2, O3, O4, P, Pd, Pdd, L, Ld, Ldd), Hddd_exact, "numpy")
fhd = sp.lambdify((a, H, O, O1, O2, P, L), Hdot_code, "numpy")
fhdd = sp.lambdify((a, H, Hd, O, O1, O2, O3, P, Pd, L, Ld), Hdotdot_code, "numpy")


def sd(lna, y):
    return CubicSpline(lna, y)(lna, 1)


print("\nnumeric (runs/A0 npz), errors = max over 1e-3<=a<=0.95 of |X - dHdotdot/dtau_spline| / max(|FD|, adotoa^4)")
print(f"{'model':24s} {'stored':>9s} {'code(eval)':>10s} {'code+fdLO':>10s} {'exact,0s':>9s} {'exact+fd':>9s}"
      f" {'ctrl Hdd':>9s}  max|Ldd_fd|/H^4 max|a^4O4_fd|")
models = ["5_quint_exponential", "5_quint_cubic", "5_quint_axion", "5_quint_inverse_power", "5_quint_quadratic",
          "5_ADE_1", "5_Kmouflage_1", "5_Kmouflage_2", "5_Kmouflage_3", "5_Kmimic_1", "5_Kmimic_2", "5_Kmimic_3",
          "5_scg"]
rows = {}
for m in models:
    d = np.load(RUNS / f"{m}.npz", allow_pickle=True)
    A = d["eft_a_grid"]
    lna = np.log(A)
    g = {k[4:]: d[k] for k in d.files if k.startswith("eft_")}
    Hh = g["adotoa"]
    Ldd_fd = Hh * sd(lna, g["EFTLambdadot"]) - 2 * Hh * g["EFTLambdadot"]
    O4_fd = sd(lna, g["EFTOmegaPPP"]) / A
    args0 = (A, Hh, g["Hdot"], g["Hdotdot"], g["EFTOmegaV"], g["EFTOmegaP"], g["EFTOmegaPP"], g["EFTOmegaPPP"],
             g["EFTOmegaPPPP"], g["gpresm_t"], g["gpresdotm_t"], g["gpresdotdotm_t"], g["EFTLambda"],
             g["EFTLambdadot"], g["EFTLambdadotdot"])
    argsf = list(args0)
    argsf[8] = O4_fd
    argsf[14] = Ldd_fd
    fd = Hh * sd(lna, g["Hdotdot"])
    sel = (A >= 1e-3) & (A <= 0.95)
    scale = np.maximum(np.abs(fd), Hh**4)[sel]

    def err(x):
        return float(np.max(np.abs(np.asarray(x)[sel] - fd[sel]) / scale))

    # controls: Hdot, Hdotdot recomputed from the coded formulas equal the stored ones
    hd_re = fhd(A, Hh, g["EFTOmegaV"], g["EFTOmegaP"], g["EFTOmegaPP"], g["gpresm_t"], g["EFTLambda"])
    hdd_re = fhdd(A, Hh, g["Hdot"], g["EFTOmegaV"], g["EFTOmegaP"], g["EFTOmegaPP"], g["EFTOmegaPPP"],
                  g["gpresm_t"], g["gpresdotm_t"], g["EFTLambda"], g["EFTLambdadot"])
    fd_hdd = Hh * sd(lna, g["Hdot"])
    ctrl_hdd = float(np.max(np.abs(g["Hdotdot"][sel] - fd_hdd[sel]) / np.maximum(np.abs(fd_hdd), Hh**3)[sel]))
    r = {
        "stored": err(g["Hdotdotdot"]),
        "code_eval_vs_stored": float(np.max(np.abs(fx(*args0) - g["Hdotdotdot"])[sel] / scale)),
        "code_with_fd_Ldd_O4": err(fx(*argsf)),
        "exact_with_stored_zeros": err(fe(*args0)),
        "exact_with_fd": err(fe(*argsf)),
        "ctrl_Hdotdot_vs_fd": ctrl_hdd,
        "recompute_Hdot_vs_stored": float(np.max(np.abs(hd_re - g["Hdot"])[sel] / (Hh**2)[sel])),
        "recompute_Hdotdot_vs_stored": float(np.max(np.abs(hdd_re - g["Hdotdot"])[sel] / (Hh**3)[sel])),
        "max|Ldd_fd|/H^4": float(np.max(np.abs(Ldd_fd[sel]) / (Hh**4)[sel])),
        "max|a^4 O4_fd|": float(np.max(np.abs(A**4 * O4_fd)[sel])),
    }
    rows[m] = r
    print(f"{m:24s} {r['stored']:9.2e} {r['code_eval_vs_stored']:10.1e} {r['code_with_fd_Ldd_O4']:10.2e}"
          f" {r['exact_with_stored_zeros']:9.2e} {r['exact_with_fd']:9.2e} {r['ctrl_Hdotdot_vs_fd']:9.1e}"
          f"  {r['max|Ldd_fd|/H^4']:9.2e} {r['max|a^4 O4_fd|']:9.2e}   (recompute Hdot/Hdotdot vs stored:"
          f" {r['recompute_Hdot_vs_stored']:.0e}/{r['recompute_Hdotdot_vs_stored']:.0e})")
out = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit/work/verify/hderivs/hddd_sympy.json")
out.write_text(json.dumps(rows, indent=1))

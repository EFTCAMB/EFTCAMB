"""Verify (mass-stability-verdict-unreliable): independent numerical evaluation of EFT_mu1/EFT_mu2 for a model.

Inputs come from the reference run npz (get_eft_functions cache on a = logspace(-8, 0, 401), no mass computation in
that path).  The Fortran expressions of EFTCAMBModelComputeMassFactors are transcribed automatically (same parser as
massverd_formula.py) and evaluated in float64 and in mpmath (50 digits) under these input variants:
  V0  code-like: EFTcdotdot = -2 adotoa EFTcdot (what compute_additional_derivs returns in A0, F-009)
  V1  EFTcdotdot from a spline derivative of EFTcdot (correct second derivative), EFTcPP as coded
  V2  V1 + EFTcPP = d^2 EFTc/da^2 (dimensionally consistent expression)
  V3  V2 + Hdot, Hdotdot replaced by their dust+DE-only values (dH/dt = -rho_d/2 - c, d2H/dt2 = 3 H rho_d/2 - dc/dt),
      i.e. radiation and neutrinos removed from the background derivatives seen by the formula
  V4  V1 + dust-only Hdot, Hdotdot (EFTcPP as coded)
Only models with Omega = gamma_i = 0 (minimally coupled) are meaningful for V3/V4.

Run: python massverd_numeric.py [npz] [amin]
"""

import importlib.util
import json
import re
import sys

import mpmath as mp
import numpy as np
from scipy.interpolate import CubicSpline

RUNS = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/runs/A0"
NPZ = sys.argv[1] if len(sys.argv) > 1 else f"{RUNS}/5_quint_exponential.npz"
AMIN = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-2
spec = importlib.util.spec_from_file_location("mf", __file__.replace("massverd_numeric.py", "massverd_formula_lib.py"))
mf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mf)

mp.mp.dps = 50
d = np.load(NPZ)
meta = json.loads(str(d["meta"]))
a = d["eft_a_grid"]
g = {k[4:]: d[k] for k in d.files if k.startswith("eft_")}
H = g["adotoa"]
# correct EFTcdotdot: dEFTcdot/dtau = adotoa * dEFTcdot/dln a; EFTcdotdot = c'' a^2 = dEFTcdot/dtau - 2 adotoa EFTcdot
la = np.log(a)
sp_cd = CubicSpline(la, g["EFTcdot"])
cdd_true = H * sp_cd(la, 1) - 2 * H * g["EFTcdot"]
cdd_code = -2 * H * g["EFTcdot"]

stmts = mf.transcribe(mf.DEFAULT_SRC, numeric=True)
stmts_fix = [("EFTcPP = (EFTcdotdot + (3*adotoa - Hdot/adotoa)*EFTcdot + 2*adotoa**2*EFTc)/(a**2*adotoa**2)"
              if s.startswith("EFTcPP") else s) for s in stmts]

# dust+DE only background derivatives (cosmic time): dH/dt = -rho_d/2 - c, d2H/dt2 = 3 H rho_d/2 - dc/dt
Hc = H / a                                   # physical H
rho_d = (g["grhob_t"] + g["grhoc_t"]) / a**2
cphys = g["EFTc"] / a**2
ct = g["EFTcdot"] / a**3                      # dc/dt = (dc/dtau)/a, dc/dtau = EFTcdot/a^2
Hd_dust = -rho_d / 2 - cphys
Hdd_dust = 1.5 * Hc * rho_d - ct
Hdot_dust = a**2 * (Hc**2 + Hd_dust)
Hdotdot_dust = a**3 * (2 * Hc**3 + 4 * Hc * Hd_dust + Hdd_dust)
Hd_code = (g["Hdot"] - H**2) / a**2
print("relative size of the non-dust part of dH/dt, (dH/dt - dust value)/H^2 at a=0.01,0.1,0.5,1:",
      [float(((Hd_code - Hd_dust) / Hc**2)[np.searchsorted(a, x)]) for x in (0.01, 0.1, 0.5, 1.0)])
print("c/H^2 at a=0.01,0.1,0.5,1:", [float((cphys / Hc**2)[np.searchsorted(a, x)]) for x in (0.01, 0.1, 0.5, 1.0)])


def evaluate(statements, overrides, use_mp=False):
    sel = np.where(a > AMIN)[0]
    out = {"mu1": np.full(len(a), np.nan), "mu2": np.full(len(a), np.nan)}
    for i in sel:
        env = {"sqrt": (mp.sqrt if use_mp else np.sqrt)}
        for k, v in g.items():
            if v.ndim == 1:
                env[k] = mp.mpf(float(v[i])) if use_mp else float(v[i])
        env["h0_Mpc"] = mp.mpf(float(H[-1])) if use_mp else float(H[-1])   # adotoa at a = 1 = H0 [1/Mpc]
        for k, v in overrides.items():
            env[k] = mp.mpf(float(v[i])) if use_mp else float(v[i])
        env["a"] = mp.mpf(float(a[i])) if use_mp else float(a[i])
        for s in statements:
            lhs, rhs = s.split("=", 1)
            env[lhs.strip()] = eval(rhs, env)  # noqa: S307
        out["mu1"][i] = float(env["EFT_mu1"]) * (a[i] / H[i]) ** 2
        out["mu2"][i] = float(env["EFT_mu2"]) * (a[i] / H[i]) ** 2
    return out


variants = {
    "V0 code-like cdd": (stmts, {"EFTcdotdot": cdd_code}),
    "V1 correct cdd": (stmts, {"EFTcdotdot": cdd_true}),
    "V2 V1+EFTcPP fixed": (stmts_fix, {"EFTcdotdot": cdd_true}),
    "V3 V2+dust-only Hdot,Hdotdot": (stmts_fix, {"EFTcdotdot": cdd_true, "Hdot": Hdot_dust, "Hdotdot": Hdotdot_dust}),
    "V4 V1+dust-only Hdot,Hdotdot": (stmts, {"EFTcdotdot": cdd_true, "Hdot": Hdot_dust, "Hdotdot": Hdotdot_dust}),
}
res = {}
for name, (st, ov) in variants.items():
    r = evaluate(st, ov)
    m = np.fmin(r["mu1"], r["mu2"])
    sel = a > AMIN
    imin = np.nanargmin(np.where(sel, m, np.inf))
    nan = int(np.sum(np.isnan(m[sel])))
    nrej = int(np.sum(m[sel] < -100))
    res[name] = {"min_mu_over_H2": float(np.nanmin(m[sel])), "a_at_min": float(a[imin]),
                 "max_abs_mu_over_H2": float(np.nanmax(np.abs(np.r_[r["mu1"][sel], r["mu2"][sel]]))),
                 "n_points_below_-100": nrej, "n_nan": nan, "n_points": int(sel.sum())}
    print(f"{name:32s} min mu/H^2 = {res[name]['min_mu_over_H2']:+.4g} at a = {a[imin]:.4g}; max|mu|/H^2 = "
          f"{res[name]['max_abs_mu_over_H2']:.3g}; points with mu/H^2 < -100: {nrej}/{int(sel.sum())}; NaN: {nan}")
    if name.startswith(("V1", "V3")):
        rm_ = evaluate(st, ov, use_mp=True)
        mm = np.fmin(rm_["mu1"], rm_["mu2"])
        rel = np.nanmax(np.abs(mm[sel] - m[sel]) / np.maximum(np.abs(mm[sel]), 1e-300))
        print(f"    mpmath(50 digits) vs float64 max relative difference: {rel:.2e}; mpmath min = {np.nanmin(mm[sel]):+.4g}")
        res[name]["mp_rel_diff"] = float(rel)
    for x in (0.02, 0.1, 0.3, 0.5, 1.0):
        j = np.searchsorted(a, x)
        print(f"    a={a[j]:.3g}: mu1/H^2={r['mu1'][j]:+.4g} mu2/H^2={r['mu2'][j]:+.4g}")
print("RESULT " + json.dumps(res))

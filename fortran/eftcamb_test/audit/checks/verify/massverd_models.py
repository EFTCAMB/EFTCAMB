"""Verify (mass-stability-verdict-unreliable): EFT_mu1/2 of several models evaluated independently of the Fortran.

Inputs: reference A0 run npz (get_eft_functions cache, a = logspace(-8,0,401); this path does not call
compute_additional_derivs or compute_mass_factors, so second derivatives are the models' own analytic values; for
full-mapping models EFTcdotdot is not filled there and is taken from a spline derivative of EFTcdot).
The Fortran mu1/mu2 block is transcribed automatically and evaluated in float64 and in mpmath (60 digits):
  A0like : Gamma{1,2,3}PP = 0 and EFTcdotdot = -2 adotoa EFTcdot (what compute_additional_derivs returns, F-009)
  own    : models' own second derivatives (what a fix of F-009 alone gives), EFTcPP as coded
  own+cPP: own + EFTcPP = d^2 EFTc/da^2 (dimensionally consistent expression)
Prints min over a > 0.01 of min(mu1, mu2)/H^2 (H physical) and the verdict at rate 10 (reject if < -100).

Run: python massverd_models.py model [model ...]
"""

import importlib.util
import json
import sys

import mpmath as mp
import numpy as np
from scipy.interpolate import CubicSpline

RUNS = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/runs/A0"
AMIN = 1e-2
spec = importlib.util.spec_from_file_location("mf", __file__.replace("massverd_models.py", "massverd_formula_lib.py"))
mf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mf)
mp.mp.dps = 60
STM = mf.transcribe(mf.DEFAULT_SRC, numeric=True)
FIX = "EFTcPP = (EFTcdotdot + (3*adotoa - Hdot/adotoa)*EFTcdot + 2*adotoa**2*EFTc)/(a**2*adotoa**2)"
STM_FIX = [FIX if s.startswith("EFTcPP") else s for s in STM]


def evaluate(a, g, statements, overrides, idx, use_mp):
    mu1 = np.full(len(a), np.nan)
    mu2 = np.full(len(a), np.nan)
    conv = (lambda x: mp.mpf(float(x))) if use_mp else float
    for i in idx:
        env = {"sqrt": mp.sqrt if use_mp else np.sqrt}
        for k, v in g.items():
            if v.ndim == 1:
                env[k] = conv(v[i])
        for k, v in overrides.items():
            env[k] = conv(v[i])
        env["a"] = conv(a[i])
        env["h0_Mpc"] = conv(g["adotoa"][-1])
        try:
            for s in statements:
                lhs, rhs = s.split("=", 1)
                env[lhs.strip()] = eval(rhs, env)  # noqa: S307
            mu1[i] = float(env["EFT_mu1"]) * (a[i] / g["adotoa"][i]) ** 2
            mu2[i] = float(env["EFT_mu2"]) * (a[i] / g["adotoa"][i]) ** 2
        except (ZeroDivisionError, ValueError):
            pass
    return mu1, mu2


def run(model):
    d = np.load(f"{RUNS}/{model}.npz")
    a = d["eft_a_grid"]
    g = {k[4:]: d[k].astype(float) for k in d.files if k.startswith("eft_") and d[k].ndim == 1}
    with np.errstate(all="ignore"):
        H = g["adotoa"]
        cdd = g["EFTcdotdot"].copy()
        if not np.any(cdd != 0):
            cdd = H * CubicSpline(np.log(a), g["EFTcdot"])(np.log(a), 1) - 2 * H * g["EFTcdot"]
        idx = np.where(a > AMIN)[0]
        variants = {
            "A0like": (STM, {"EFTGamma1PP": 0 * a, "EFTGamma2PP": 0 * a, "EFTGamma3PP": 0 * a,
                             "EFTcdotdot": -2 * H * g["EFTcdot"]}),
            "own": (STM, {"EFTcdotdot": cdd}),
            "own+cPP": (STM_FIX, {"EFTcdotdot": cdd}),
        }
        out = {"model": model, "c_over_H2_at_0.01_0.1_1": [float((g["EFTc"] / H**2)[np.searchsorted(a, x)])
                                                            for x in (0.01, 0.1, 1.0)]}
        for name, (st, ov) in variants.items():
            row = {}
            for prec, use_mp in (("f64", False), ("mp60", True)):
                m1, m2 = evaluate(a, g, st, ov, idx, use_mp)
                m = np.fmin(m1, m2)[idx]
                nn = int(np.sum(np.isnan(np.fmin(m1, m2)[idx]) | np.isnan(m1[idx]) | np.isnan(m2[idx])))
                mn = float(np.nanmin(m)) if np.any(np.isfinite(m)) else float("nan")
                row[prec] = {"min": mn, "a_min": float(a[idx][np.nanargmin(m)]) if np.isfinite(mn) else None,
                             "nan": nn, "reject": bool(mn < -100 or nn > 0)}
            out[name] = row
    return out


if __name__ == "__main__":
    models = sys.argv[1:] or ["5_quint_exponential"]
    allres = []
    for mdl in models:
        r = run(mdl)
        allres.append(r)
        print(f"{mdl}: c/H^2 at a=0.01,0.1,1 = {', '.join(f'{x:.2g}' for x in r['c_over_H2_at_0.01_0.1_1'])}")
        for name in ("A0like", "own", "own+cPP"):
            f, m = r[name]["f64"], r[name]["mp60"]
            print(f"    {name:8s} float64: min mu/H^2 = {f['min']:+.4g} (a={f['a_min']}) NaN={f['nan']} "
                  f"{'REJECT' if f['reject'] else 'accept'} | mpmath60: min = {m['min']:+.4g} (a={m['a_min']}) "
                  f"NaN={m['nan']} {'REJECT' if m['reject'] else 'accept'}")
    print("RESULT " + json.dumps(allres))

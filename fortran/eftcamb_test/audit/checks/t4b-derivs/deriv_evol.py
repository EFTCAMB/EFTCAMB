"""Track t4b-derivs: checks that need get_scale_evolution (fields only filled inside derivs).

1) background consistency: every background field of the cache filled in derivs (ev_*) against get_eft_functions at
   the same scale factor (evst_*);
2) derivative relations on the evolution stencil (non-uniform 9-point FD in ln a, d/dtau = adotoa d/dln a):
   background-only Einstein factors (EFTeomXdot, EFTeomVdot) and perturbation fields (EFTeomNdot, pidot, ...).

Usage: python deriv_evol.py DIR [--json OUT]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

BG_FIELDS = [
    "grhom_t", "gpresm_t", "gpresdotm_t", "gpresdotdotm_t", "grhonu_tot", "gpinu_tot", "grhonudot_tot",
    "gpinudot_tot", "gpinudotdot_tot", "adotoa", "Hdot", "Hdotdot", "Hdotdotdot", "EFTOmegaV", "EFTOmegaP",
    "EFTOmegaPP", "EFTOmegaPPP", "EFTOmegaPPPP", "EFTc", "EFTcdot", "EFTcdotdot", "EFTLambda", "EFTLambdadot",
    "EFTLambdadotdot", "EFTGamma1V", "EFTGamma1P", "EFTGamma1PP", "EFTGamma2V", "EFTGamma2P", "EFTGamma2PP",
    "EFTGamma2PPP", "EFTGamma3V", "EFTGamma3P", "EFTGamma3PP", "EFTGamma3PPP", "EFTGamma3PPPP", "EFTGamma4V",
    "EFTGamma4P", "EFTGamma4PP", "EFTGamma5V", "EFTGamma5P", "EFTGamma6V", "EFTGamma6P", "alphaB", "alphaBdot",
    "alphaT", "alphaTdot", "alphaK", "alphaKdot", "alphaM", "alphaMdot", "Meff2", "phi_scf", "dphi", "ddphi",
    "dddphi", "grhoq", "gpresq", "grhodotq", "gpresdotq", "grhov_t", "gpiv_t",
]

# (stored, base, kind, order, is_background_only)
EV_RELATIONS = [
    ("EFTeomXdot", "EFTeomX", "tau", 1, True),
    ("EFTeomVdot", "EFTeomV", "tau", 1, True),
    ("EFTeomNdot", "EFTeomN", "tau", 1, False),
    ("sijdot", "sij", "tau", 1, False),
    ("pidot", "pi", "tau", 1, False),
    ("pidotdot", "pidot", "tau", 2, False),
    ("delta_phidot", "delta_phi", "tau", 1, False),
    ("delta_phidotdot", "delta_phidot", "tau", 2, False),
    ("PsiDot", "Psi", "tau", 1, False),
    ("PhiDot", "Phi", "tau", 1, False),
    ("sigmadot", "sigma", "tau", 1, False),
    ("sigmadotdot", "sigmadot", "tau", 2, False),
    ("etakdot", "etak", "tau", 1, False),
    ("dz", "z", "tau", 1, False),
]


def fd_weights(x: np.ndarray, x0: float, m: int = 1) -> np.ndarray:
    """Fornberg finite-difference weights for the m-th derivative at x0 on nodes x."""
    n = len(x)
    c = np.zeros((n, m + 1))
    c1, c4 = 1.0, x[0] - x0
    c[0, 0] = 1.0
    for i in range(1, n):
        mn = min(i, m)
        c2, c5, c4 = 1.0, c4, x[i] - x0
        for j in range(i):
            c3 = x[i] - x[j]
            c2 *= c3
            if j == i - 1:
                for k in range(mn, 0, -1):
                    c[i, k] = c1 * (k * c[i - 1, k - 1] - c5 * c[i - 1, k]) / c2
                c[i, 0] = -c1 * c5 * c[i - 1, 0] / c2
            for k in range(mn, 0, -1):
                c[j, k] = (c4 * c[j, k] - k * c[j, k - 1]) / c3
            c[j, 0] = c4 * c[j, 0] / c3
        c1 = c2
    return c[:, m]


def deriv_lna(xs: np.ndarray, f: np.ndarray, npts: int = 9) -> np.ndarray:
    """xs, f: (n_c, 9). d f/d ln a at the stencil centre using npts central nodes."""
    lo = (9 - npts) // 2
    out = np.empty(xs.shape[0])
    for i in range(xs.shape[0]):
        x = xs[i, lo : lo + npts]
        w = fd_weights(x, xs[i, 4])
        out[i] = w @ f[i, lo : lo + npts]
    return out


def load(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as f:
        d = {k: f[k] for k in f.files}
    d["meta"] = json.loads(str(d["meta"]))
    return d


def analyze(path: Path) -> dict:
    d = load(path)
    res = {"model": path.stem, "background": {}, "relations": {}}
    if "ev_a" not in d:
        return res
    t_on = max(d["meta"]["pert_turn_on_results"], d["meta"]["back_turn_on"])
    for ik, k in enumerate(d["ev_k"]):
        a = d["ev_a"][ik]
        valid = np.all(a > 0, axis=1)
        h = float(d["ev_h"])
        on = valid & (a[:, 0] > t_on * np.exp(h))
        # 1) background: derivs vs getter
        for f in BG_FIELDS:
            x, y = d[f"ev_{f}"][ik][:, 4], d[f"evst_{f}"][ik][:, 4]
            if not (np.any(x != 0) or np.any(y != 0)):
                continue
            sc = np.maximum(np.abs(x), np.abs(y))
            sc = np.where(sc > 0, sc, 1.0)
            e = np.where(on, np.abs(x - y) / sc, 0.0)
            i = int(np.argmax(e))
            prev = res["background"].get(f, {"max_rel": -1})
            if e[i] > prev["max_rel"]:
                res["background"][f] = {"max_rel": float(e[i]), "a": float(a[i, 4]), "k": float(k),
                                        "derivs": float(x[i]), "getter": float(y[i])}
        # 2) derivative relations along the trajectory
        xs = np.log(np.where(a > 0, a, 1.0))
        H = d["ev_adotoa"][ik][:, 4]
        for stored, base, kind, order, bg in EV_RELATIONS:
            S = d[f"ev_{stored}"][ik][:, 4]
            X = d[f"ev_{base}"][ik]
            if not (np.any(S[on] != 0) or np.any(X[on] != 0)):
                continue
            D9 = deriv_lna(xs[on], X[on], 9)
            D5 = deriv_lna(xs[on], X[on], 5)
            pred9, pred5 = H[on] * D9, H[on] * D5
            sc = np.maximum.reduce([np.abs(S[on]), np.abs(pred9), H[on] * np.abs(X[on, 4])])
            n = len(sc)
            sc = np.array([np.max(sc[max(0, i - 2) : i + 3]) for i in range(n)])
            sc = np.where(sc > 0, sc, 1.0)
            err = np.abs(S[on] - pred9) / sc
            spread = np.abs(pred9 - pred5) / sc
            excess = np.maximum(err - 3 * spread, 0)
            i = int(np.argmax(excess))
            key = f"{stored}@k={k:g}"
            tol = 1e-6 if order == 1 else 1e-4
            res["relations"][key] = {
                "max_excess": float(excess[i]), "a": float(a[on][i, 4]), "err": float(err[i]),
                "spread": float(spread[i]), "median_err": float(np.median(err)), "median_spread": float(np.median(spread)),
                "background_only": bg, "tol": tol, "pass": bool(excess.max() < tol),
            }
    return res


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("dir")
    p.add_argument("--json")
    p.add_argument("--model", default="")
    a = p.parse_args()
    out = []
    for f in sorted(Path(a.dir).glob(f"{a.model}*.npz")):
        r = analyze(f)
        out.append(r)
        bad_bg = {k: v for k, v in r["background"].items() if v["max_rel"] > 1e-10}
        print(f"== {r['model']}")
        for k, v in bad_bg.items():
            print(f"   BG-MISMATCH {k:16s} rel={v['max_rel']:.2e} @a={v['a']:.3e} derivs={v['derivs']:.6e} "
                  f"getter={v['getter']:.6e}")
        for k, v in r["relations"].items():
            flag = "PASS" if v["pass"] else "FAIL"
            print(f"   {flag} {k:22s} excess={v['max_excess']:.2e} @a={v['a']:.2e} err={v['err']:.2e} "
                  f"spread={v['spread']:.2e} med_err={v['median_err']:.1e} med_spread={v['median_spread']:.1e}"
                  f"{' [bg]' if v['background_only'] else ''}")
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

"""Track t4b-derivs: check derivative relations of the EFTCAMB timestep cache from deriv_worker.py stencil output.

Conventions (from 03_abstract_EFTCAMB_cache.f90 and the model code):
  ' = d/dtau (conformal), d/dtau = adotoa d/dln a;  P suffix = d/da.
  kind "a"     : stored = dX/da                        = D(X)/a
  kind "tau"   : stored = dX/dtau                      = adotoa D(X)
  kind "tau_a2": stored = a^2 d(X/a^2)/dtau            = adotoa (D(X) - 2 X)   (X is an a^2-weighted quantity)
with D = d/dln a evaluated by 8th-order central differences at three step sizes.

Usage: python deriv_analyze.py DIR [--json OUT] [--model NAME] [--verbose]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

W8 = np.array([1 / 280, -4 / 105, 1 / 5, -4 / 5, 0.0, 4 / 5, -1 / 5, 4 / 105, -1 / 280])
TOL1, TOLH = 1e-6, 1e-4

# (stored field, base field, kind, order of the stored derivative)
RELATIONS = [
    ("Hdot", "adotoa", "tau", 1),
    ("Hdotdot", "Hdot", "tau", 2),
    ("Hdotdotdot", "Hdotdot", "tau", 3),
    ("EFTOmegaP", "EFTOmegaV", "a", 1),
    ("EFTOmegaPP", "EFTOmegaP", "a", 2),
    ("EFTOmegaPPP", "EFTOmegaPP", "a", 3),
    ("EFTOmegaPPPP", "EFTOmegaPPP", "a", 4),
    ("EFTcdot", "EFTc", "tau_a2", 1),
    ("EFTcdotdot", "EFTcdot", "tau_a2", 2),
    ("EFTLambdadot", "EFTLambda", "tau_a2", 1),
    ("EFTLambdadotdot", "EFTLambdadot", "tau_a2", 2),
    ("EFTGamma1P", "EFTGamma1V", "a", 1),
    ("EFTGamma1PP", "EFTGamma1P", "a", 2),
    ("EFTGamma2P", "EFTGamma2V", "a", 1),
    ("EFTGamma2PP", "EFTGamma2P", "a", 2),
    ("EFTGamma2PPP", "EFTGamma2PP", "a", 3),
    ("EFTGamma3P", "EFTGamma3V", "a", 1),
    ("EFTGamma3PP", "EFTGamma3P", "a", 2),
    ("EFTGamma3PPP", "EFTGamma3PP", "a", 3),
    ("EFTGamma3PPPP", "EFTGamma3PPP", "a", 4),
    ("EFTGamma4P", "EFTGamma4V", "a", 1),
    ("EFTGamma4PP", "EFTGamma4P", "a", 2),
    ("EFTGamma5P", "EFTGamma5V", "a", 1),
    ("EFTGamma6P", "EFTGamma6V", "a", 1),
    ("alphaBdot", "alphaB", "tau", 1),
    ("alphaTdot", "alphaT", "tau", 1),
    ("alphaKdot", "alphaK", "tau", 1),
    ("alphaMdot", "alphaM", "tau", 1),
    ("gpresdotm_t", "gpresm_t", "tau_a2", 1),
    ("gpresdotdotm_t", "gpresdotm_t", "tau_a2", 2),
    ("grhonudot_tot", "grhonu_tot", "tau_a2", 1),
    ("gpinudot_tot", "gpinu_tot", "tau_a2", 1),
    ("gpinudotdot_tot", "gpinudot_tot", "tau_a2", 2),
    ("dphi", "phi_scf", "tau", 1),
    ("ddphi", "dphi", "tau", 2),
    ("dddphi", "ddphi", "tau", 3),
    ("grhodotq", "grhoq", "tau_a2", 1),
    ("gpresdotq", "gpresq", "tau_a2", 1),
    ("sijdot", "sij", "tau", 1),
]


# physical floor of the stored derivative: phys0 * rate**order, rate = 1/a ("a" kinds) or adotoa ("tau" kinds)
# phys0: 1 for dimensionless functions (Omega, Gamma_i, alpha_i, phi in Planck units), adotoa for adotoa,
# adotoa**2 for a^2-weighted densities/pressures and the EFT c, Lambda.
FAMILY = {
    "adotoa": "H", "Hdot": "H", "Hdotdot": "H",
    "EFTc": "rho", "EFTcdot": "rho", "EFTLambda": "rho", "EFTLambdadot": "rho",
    "gpresm_t": "rho", "gpresdotm_t": "rho", "grhonu_tot": "rho", "gpinu_tot": "rho", "gpinudot_tot": "rho",
    "grhoq": "rho", "gpresq": "rho", "sij": None,
}


def phys_floor(d: dict, base: str, kind: str, order: int) -> np.ndarray:
    a = np.exp(d["xc"])
    H = d["st_adotoa"][0, 4]
    fam = FAMILY.get(base, "dimless")
    if fam is None:
        return np.zeros_like(a)
    phys0 = {"dimless": 1.0, "H": H, "rho": H**2}[fam]
    # order of the stored derivative relative to the undifferentiated family member
    n = order if fam == "dimless" else {"H": 1, "rho": 1}[fam] * order
    rate = (1.0 / a) if kind == "a" else H
    return phys0 * rate**n


def load(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as f:
        d = {k: f[k] for k in f.files}
    d["meta"] = json.loads(str(d["meta"]))
    return d


def dlna(st: np.ndarray, hs: np.ndarray) -> np.ndarray:
    """st: (n_h, 9, n_c) -> (n_h, n_c) d/dln a."""
    return np.einsum("j,hjc->hc", W8, st) / hs[:, None]


def window_max(x: np.ndarray, w: int = 3) -> np.ndarray:
    n = len(x)
    return np.array([np.max(x[max(0, i - w) : min(n, i + w + 1)]) for i in range(n)])


def predict(d: dict, base: str, kind: str):
    hs = d["hs"]
    X = d[f"st_{base}"]
    a = np.exp(d["xc"])
    H = d["st_adotoa"][0, 4]
    D = dlna(X, hs)  # (n_h, n_c)
    Xc = X[0, 4]
    if kind == "a":
        pred = D / a
        natural = np.abs(Xc) / a
    elif kind == "tau":
        pred = H * D
        natural = H * np.abs(Xc)
    elif kind == "tau_a2":
        pred = H * (D - 2.0 * Xc)
        natural = H * np.abs(Xc)
    else:
        raise ValueError(kind)
    return pred, natural


def check_relation(d: dict, stored: str, base: str, kind: str, order: int, mask: np.ndarray | None = None) -> dict | None:
    if f"st_{stored}" not in d or f"st_{base}" not in d:
        return None
    S = d[f"st_{stored}"][0, 4]
    X = d[f"st_{base}"]
    if not np.any(X != 0) and not np.any(S != 0):
        return None  # field not used by this model
    pred, natural = predict(d, base, kind)
    pm = pred[1]
    scale = window_max(np.maximum.reduce([np.abs(S), np.abs(pm), natural]))
    scale = np.where(scale > 0, scale, 1.0)
    err = np.abs(S - pm) / scale
    spread = np.max(np.abs(pred - pm[None, :]), axis=0) / scale
    good = np.isfinite(err) & np.isfinite(spread)
    if mask is not None:
        good &= mask
    tol = TOL1 if order == 1 else TOLH
    resolved = good & (spread < 0.1 * tol)
    excess = np.where(good, np.maximum(err - 3 * spread, 0.0), 0.0)
    floor = phys_floor(d, base, kind, order)
    sc_phys = np.maximum(scale, floor)
    err_p = np.abs(S - pm) / sc_phys
    spread_p = np.max(np.abs(pred - pm[None, :]), axis=0) / sc_phys
    excess_p = np.where(good, np.maximum(err_p - 3 * spread_p, 0.0), 0.0)
    i_p = int(np.argmax(excess_p))
    a = np.exp(d["xc"])
    i_err = int(np.argmax(np.where(resolved, err, -1))) if resolved.any() else -1
    i_exc = int(np.argmax(excess))
    # alternative hypotheses for tau-kind fields (to establish conventions)
    alt = {}
    if kind in ("tau", "tau_a2"):
        other = "tau_a2" if kind == "tau" else "tau"
        p2, n2 = predict(d, base, other)
        sc2 = window_max(np.maximum.reduce([np.abs(S), np.abs(p2[1]), n2]))
        sc2 = np.where(sc2 > 0, sc2, 1.0)
        e2 = np.abs(S - p2[1]) / sc2
        alt[other] = float(np.nanmax(np.where(good, e2, 0.0)))
    return {
        "stored": stored,
        "base": base,
        "kind": kind,
        "order": order,
        "tol": tol,
        "n_points": int(good.sum()),
        "n_resolved": int(resolved.sum()),
        "max_err_resolved": float(err[i_err]) if i_err >= 0 else None,
        "a_max_err_resolved": float(a[i_err]) if i_err >= 0 else None,
        "max_excess": float(excess[i_exc]),
        "a_max_excess": float(a[i_exc]),
        "spread_at_max_excess": float(spread[i_exc]),
        "err_at_max_excess": float(err[i_exc]),
        "median_err": float(np.median(err[good])) if good.any() else None,
        "max_excess_phys": float(excess_p[i_p]),
        "a_max_excess_phys": float(a[i_p]),
        "pass_local": bool((excess[good] < tol).all()),
        "pass": bool((excess_p[good] < tol).all()),
        "alt": alt,
    }


def generic_relations(d: dict, mask: np.ndarray | None = None) -> dict:
    """Background identities that must hold for any EFT model (flat universe)."""
    g = lambda f: d[f"st_{f}"][0, 4]  # noqa: E731
    a = np.exp(d["xc"])
    H, Hd = g("adotoa"), g("Hdot")
    Om, OmP, OmPP = g("EFTOmegaV"), g("EFTOmegaP"), g("EFTOmegaPP")
    c, L = g("EFTc"), g("EFTLambda")
    rho, p = g("grhom_t"), g("gpresm_t")
    out = {}

    def rel(name, lhs, rhs, extra_scale=None):
        sc = np.maximum(np.abs(lhs), np.abs(rhs))
        if extra_scale is not None:
            sc = np.maximum(sc, extra_scale)
        sc = np.where(sc > 0, sc, 1.0)
        e = np.abs(lhs - rhs) / sc
        m = np.isfinite(e) if mask is None else (np.isfinite(e) & mask)
        i = int(np.argmax(np.where(m, e, -1)))
        out[name] = {"max_err": float(e[i]), "a_at_max": float(a[i]), "median": float(np.median(e[m])), "pass": bool(e[i] < TOL1)}

    # first Friedmann equation in EFT form: (1+Omega+a Omega') H^2 = (rho_m + 2c - Lambda)/3
    rel("Friedmann_EFT", (1 + Om + a * OmP) * H**2, (rho + 2 * c - L) / 3, extra_scale=H**2)
    # second Friedmann equation in EFT form
    lhs2 = (1 + Om) * (2 * Hd + H**2) + a * OmP * (Hd + 2 * H**2) + a**2 * OmPP * H**2
    rel("Hdot_EFT", lhs2, -(p + L), extra_scale=np.maximum(H**2, np.abs(Hd)))
    # the same through the effective Q fluid (compute_rhoQPQ)
    rel("Friedmann_Q", (1 + Om) * H**2, (rho + g("grhoq")) / 3, extra_scale=H**2)
    rel("Hdot_Q", (1 + Om) * (2 * Hd + H**2), -(p + g("gpresq")), extra_scale=np.maximum(H**2, np.abs(Hd)))
    # designer-type background via grhov_t / gpiv_t (only where the model fills them)
    if np.any(g("grhov_t") != 0):
        rel("Friedmann_grhov", H**2, (rho + g("grhov_t")) / 3, extra_scale=H**2)
        if np.any(g("gpiv_t") != 0):
            rel("Hdot_gpiv", Hd, -0.5 * (H**2 + p + g("gpiv_t")), extra_scale=np.maximum(H**2, np.abs(Hd)))
            # DE conservation: d(grhov_t)/dtau = -adotoa (grhov_t + 3 gpiv_t)
            D = dlna(d["st_grhov_t"], d["hs"])[1]
            rel("DE_conservation", H * D, -H * (g("grhov_t") + 3 * g("gpiv_t")), extra_scale=H * np.abs(g("grhov_t")))
    # matter conservation for the a^2-weighted densities
    for f, w in (("grhob_t", 0.0), ("grhoc_t", 0.0), ("grhor_t", 1 / 3), ("grhog_t", 1 / 3)):
        D = dlna(d[f"st_{f}"], d["hs"])[1]
        rel(f"cons_{f}", H * D, -H * (1 + 3 * w) * g(f), extra_scale=H * np.abs(g(f)))
    D = dlna(d["st_grhom_t"], d["hs"])[1]
    rel("cons_grhom_t", H * D, -H * (rho + 3 * p), extra_scale=H * rho)
    if np.any(g("grhonu_tot") != 0):
        rel("cons_grhonu", g("grhonudot_tot"), -3 * H * (g("grhonu_tot") + g("gpinu_tot")), extra_scale=H * g("grhonu_tot"))
    rel("grhoa2_def", g("grhoa2"), rho / a**4)
    # adotoa against CAMB's own background (dtauda): adotoa = a H(z)/c
    Hc = d["camb_H"] / float(d["camb_c_km_s"]) * a
    rel("adotoa_vs_CAMB_H", H, Hc)
    # tau consistency: d tau / d ln a = 1/adotoa
    Dt = dlna(d["st_tau"], d["hs"])[1]
    rel("dtau_dlna", Dt * H, np.ones_like(H))
    return out


def analyze(path: Path) -> dict:
    d = load(path)
    res = {"model": path.stem, "meta": {k: v for k, v in d["meta"].items() if k != "fields"}, "relations": {},
           "generic": {}}
    a = np.exp(d["xc"])
    a_on = max(d["meta"].get("pert_turn_on_results", 0.0), d["meta"].get("back_turn_on", 0.0))
    active = a * np.exp(-4 * np.max(d["hs"])) > a_on
    res["a_on"] = a_on
    for stored, base, kind, order in RELATIONS:
        r = check_relation(d, stored, base, kind, order)
        if r is not None:
            ra = check_relation(d, stored, base, kind, order, mask=active) if active.any() else None
            r["active"] = None if ra is None else {k: ra[k] for k in ("max_excess", "a_max_excess", "max_excess_phys",
                                                                    "a_max_excess_phys", "pass", "pass_local")}
            res["relations"][stored] = r
    res["generic"] = generic_relations(d)
    res["generic_active"] = generic_relations(d, mask=active) if active.any() else {}
    return res


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("dir")
    p.add_argument("--json")
    p.add_argument("--model", default="")
    p.add_argument("-v", "--verbose", action="store_true")
    a = p.parse_args()
    out = []
    for f in sorted(Path(a.dir).glob(f"{a.model}*.npz")):
        try:
            r = analyze(f)
        except Exception as e:  # noqa: BLE001
            print(f"{f.stem}: ERROR {e!r}", file=sys.stderr)
            continue
        out.append(r)
        print(f"== {r['model']}  designer={r['meta'].get('designer')}  RGR={r['meta'].get('pert_turn_on_results')}")
        for k, v in r["relations"].items():
            flag = "PASS" if v["pass"] else "FAIL"
            if a.verbose or not v["pass"]:
                mer = v["max_err_resolved"]
                print(f"   {flag} {k:16s} phys={v['max_excess_phys']:.2e} @a={v['a_max_excess_phys']:.2e} "
                      f"local={v['max_excess']:.2e} @a={v['a_max_excess']:.2e} "
                      f"(err={v['err_at_max_excess']:.2e} spread={v['spread_at_max_excess']:.2e}) "
                      f"maxerr_res={mer if mer is None else f'{mer:.2e}'} nres={v['n_resolved']}/{v['n_points']} "
                      f"alt={ {k2: f'{v2:.1e}' for k2, v2 in v['alt'].items()} }")
        for k, v in r["generic"].items():
            flag = "PASS" if v["pass"] else "FAIL"
            if a.verbose or not v["pass"]:
                print(f"   {flag} {k:20s} max={v['max_err']:.2e} @a={v['a_at_max']:.2e} median={v['median']:.2e}")
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

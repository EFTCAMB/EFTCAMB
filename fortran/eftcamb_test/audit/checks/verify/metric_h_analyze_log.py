"""Analyse the derivs logs written by the instrumented tree (metric_h_build_instr.py / metric_h_instr_runs.py).

Record (16 float64): k, tau, a, h'(used), h'_00 (constraint), h''(ii), etak, dgrho, dphi, dphi', RSA, TC, no_phot, dgq,
etak', sigma.
For each logged k: (i) T runs: constraint violation h'_used/h'_00 - 1 along a; (ii) F run: consistency of the ii equation
with the 00 constraint, i.e. d(h'_00)/dtau (from the logged trajectory) vs h''_ii, as relative residual
r = (d h'_00/dtau - h''_ii) / max(|h''_ii|, |2 k etak|).

Usage: python metric_h_analyze_log.py <tag>   (tag e.g. 5_hdsk_TG_1)
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "metric_h"
tag = sys.argv[1] if len(sys.argv) > 1 else "5_hdsk_TG_1"
COLS = "k tau a hdot hc hdd etak dgrho dphi dphidot rsa tc nophot dgq etakdot sigma".split()
A_SAMPLES = [1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 0.1, 0.3, 1.0]
TARGETS = [1e-3, 3e-3, 1e-2, 3e-2, 0.1, 1.0, 10.0]


def load(name):
    raw = np.fromfile(OUT / f"{tag}_{name}.bin", dtype=np.float64)
    raw = raw[: raw.size // 16 * 16].reshape(-1, 16)
    return {c: raw[:, i] for i, c in enumerate(COLS)}


def per_k(d, kt):
    ks = np.unique(d["k"])
    kk = ks[np.argmin(abs(np.log(ks / kt)))]
    if abs(np.log(kk / kt)) > 0.03:
        return None, None
    m = d["k"] == kk
    o = np.argsort(d["tau"][m], kind="stable")
    return kk, {c: v[m][o] for c, v in d.items()}


def at_a(x, y, a):
    return np.interp(np.log(a), np.log(x["a"]), y, left=np.nan, right=np.nan)


res = {}
logs = {n: load(n) for n in ("F", "T", "Ticfix") if (OUT / f"{tag}_{n}.bin").exists()}
print(f"tag {tag}; a samples: " + " ".join(f"{a:g}" for a in A_SAMPLES))
for kt in TARGETS:
    out = {}
    for n, d in logs.items():
        kk, x = per_k(d, kt)
        if x is None:
            continue
        if n == "F":
            # F: h' used is h'_00 by construction; residual of the ii equation along the trajectory.
            # local quadratic fit of h'_00(tau) in windows of the sorted samples -> derivative at window centre
            t, hc, hdd = x["tau"], x["hc"], x["hdd"]
            ut, idx = np.unique(t, return_index=True)
            hc_u, hdd_u = hc[idx], hdd[idx]
            n_w = 41
            cen, dfd, hddm, scale = [], [], [], []
            for i0 in range(0, len(ut) - n_w, n_w // 2):
                sl = slice(i0, i0 + n_w)
                tc = ut[sl].mean()
                p = np.polyfit(ut[sl] - tc, hc_u[sl], 3)
                cen.append(tc)
                dfd.append(p[-2])
                hddm.append(np.interp(tc, ut[sl], hdd_u[sl]))
                scale.append(max(np.abs(hdd_u[sl]).max(), np.abs(2 * kk * x["etak"][idx][sl]).max()))
            cen, dfd, hddm, scale = map(np.array, (cen, dfd, hddm, scale))
            a_c = np.interp(cen, ut, x["a"][idx])
            r = (dfd - hddm) / scale
            samp = [float(np.interp(np.log(a), np.log(a_c), r, left=np.nan, right=np.nan)) for a in A_SAMPLES]
            out["F_ii_residual"] = samp
            out["F_sanity_hdot_eq_hc"] = float(np.nanmax(np.abs(x["hdot"] / x["hc"] - 1)))
            out["F_ii_residual_max_late"] = float(np.nanmax(np.abs(r[a_c > 1e-3]))) if np.any(a_c > 1e-3) else None
        else:
            v = x["hdot"] / x["hc"] - 1
            out[f"{n}_violation"] = [float(at_a(x, v, a)) for a in A_SAMPLES]
        out["k"] = float(kk)
    res[kt] = out
    print(f"k target {kt:g} (logged k={out.get('k', float('nan')):.5g})")
    for key, val in out.items():
        if isinstance(val, list):
            print(f"   {key:22s}: " + " ".join(f"{v:+.2e}" for v in val))
        elif key != "k":
            print(f"   {key:22s}: {val}")
(OUT / f"{tag}_loganalysis.json").write_text(json.dumps(res, indent=1))

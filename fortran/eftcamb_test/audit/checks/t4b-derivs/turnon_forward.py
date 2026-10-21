"""Track t4b-derivs: derivative relations of the timestep cache AT the EFT perturbation turn-on time a_on.

The EFT initial conditions (09_EFTCAMB_IC.f90) use the cache at a_on = EFTCAMB_pert_turn_on (after RGR).  Central
stencils around a_on straddle the switch (and, for several full-mapping models, the start of the background
table), so here the derivative of each base field is estimated with one-sided (forward) 8th-order differences in
ln a on a_on * exp(j h), j = 0..8, for h = 0.02, 0.01, 0.005, and compared with the stored derivative at a_on.
Error scale: max(|stored|, |FD|, natural scale, physical floor) as in deriv_analyze.check_relation; a relation is
'resolved' when the three step sizes agree to better than 1/3 of the error.

    python turnon_forward.py [--tree A0] [--jobs 4] [--pattern '*.ini'] OUT.json
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import harness  # noqa: E402

HS = (0.02, 0.01, 0.005)


def fwd_weights(n: int = 9) -> np.ndarray:
    """Weights w_j (j=0..n-1) with sum_j w_j f(j h) = h f'(0) + O(h^n)."""
    j = np.arange(n, dtype=float)
    V = np.vander(j, n, increasing=True).T  # V[p, j] = j**p
    rhs = np.zeros(n)
    rhs[1] = 1.0
    return np.linalg.solve(V, rhs)


def worker(cfg: dict) -> None:
    harness._exit_with_parent()
    sys.path.insert(0, cfg["tree"])
    import camb

    from deriv_analyze import FAMILY, RELATIONS

    cp = camb.read_ini(cfg["ini"], no_validate=True)
    res = camb.get_background(cp)
    a_on = float(res.Params.EFTCAMB.EFTCAMB_pert_turn_on)
    out = {"a_on": a_on, "relations": {}}
    if not (0 < a_on < 1):
        print(json.dumps(out))
        return
    W = fwd_weights()
    pts = np.unique(np.concatenate([a_on * np.exp(h * np.arange(9)) for h in HS]))
    _, arr = cp.EFTCAMB.get_eft_functions(res, pts)
    idx = {h: np.searchsorted(pts, a_on * np.exp(h * np.arange(9))) for h in HS}
    i0 = idx[HS[0]][0]
    g = {f: np.asarray(arr[f], dtype=float) for f in arr.dtype.names}
    H = g["adotoa"][i0]
    for stored, base, kind, order in RELATIONS:
        if stored not in g or base not in g:
            continue
        S, X0 = g[stored][i0], g[base][i0]
        if S == 0 and X0 == 0:
            continue
        preds = []
        for h in HS:
            D = W @ g[base][idx[h]] / h
            preds.append({"a": D / a_on, "tau": H * D, "tau_a2": H * (D - 2 * X0)}[kind])
        preds = np.array(preds)
        p = preds[-1]
        natural = abs(X0) / a_on if kind == "a" else H * abs(X0)
        fam = FAMILY.get(base, "dimless")
        rate = 1 / a_on if kind == "a" else H
        floor = 0.0 if fam is None else {"dimless": 1.0, "H": H, "rho": H**2}[fam] * rate**order
        sc = max(abs(S), abs(p), natural, floor) or 1.0
        err = abs(S - p) / sc
        spread = float(np.max(np.abs(preds - p))) / sc
        out["relations"][stored] = {"err": err, "spread": spread, "stored": S, "fd": float(p),
                                    "resolved": bool(spread < err / 3 or err < 1e-12), "order": order}
    # pointwise background identities at a_on
    Om, OmP, OmPP = g["EFTOmegaV"][i0], g["EFTOmegaP"][i0], g["EFTOmegaPP"][i0]
    c, L, rho, p_ = g["EFTc"][i0], g["EFTLambda"][i0], g["grhom_t"][i0], g["gpresm_t"][i0]
    Hd = g["Hdot"][i0]
    f1 = ((1 + Om + a_on * OmP) * H**2 - (rho + 2 * c - L) / 3) / H**2
    f2 = ((1 + Om) * (2 * Hd + H**2) + a_on * OmP * (Hd + 2 * H**2) + a_on**2 * OmPP * H**2 + p_ + L) / H**2
    Hc = float(res.hubble_parameter(1 / a_on - 1)) / (camb.constants.c / 1e3) * a_on
    out["generic"] = {"Friedmann_EFT": abs(f1), "Hdot_EFT": abs(f2), "adotoa_vs_CAMB_H": abs(H / Hc - 1)}
    print(json.dumps(out))


def run_one(tree: str, ini: Path, outdir: Path) -> dict:
    wrapper = harness._wrapper_ini(ini, outdir, ini.stem, {"EFTCAMB_feedback_level": "0", "feedback_level": "0"})
    cfg = {"tree": str(harness.tree_path(tree)), "ini": str(wrapper)}
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, __file__, "--worker", json.dumps(cfg)], cwd=outdir, env=env,
                           capture_output=True, text=True, timeout=900, preexec_fn=harness._nice)
    lines = [x for x in p.stdout.splitlines() if x.startswith("{")]
    if p.returncode != 0 or not lines:
        return {"model": ini.stem, "status": f"rc={p.returncode}", "tail": (p.stdout + p.stderr)[-800:]}
    return {"model": ini.stem, "status": "ok", **json.loads(lines[-1])}


def main() -> None:
    if len(sys.argv) > 2 and sys.argv[1] == "--worker":
        worker(json.loads(sys.argv[2]))
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--tree", default="A0")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--pattern", default="*.ini")
    a = ap.parse_args()
    outdir = Path(a.out).resolve().with_suffix("")
    outdir.mkdir(parents=True, exist_ok=True)
    inis = [i for i in harness.model_inis(a.pattern) if i.stem != "1_EFT_GR"]
    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(lambda i: run_one(a.tree, i, outdir), inis))
    Path(a.out).write_text(json.dumps(res, indent=1))
    for r in res:
        if r["status"] != "ok":
            print(f"{r['model']:30s} {r['status']}")
            continue
        bad = {k: v for k, v in r["relations"].items() if v["resolved"] and v["err"] > (1e-6 if v["order"] == 1 else 1e-4)}
        gen = {k: v for k, v in r.get("generic", {}).items() if v > 1e-6}
        print(f"{r['model']:30s} a_on={r['a_on']:.3e} n_rel={len(r['relations']):2d} "
              f"FAIL(resolved): {', '.join(f'{k}={v['err']:.1e}' for k, v in bad.items()) or '-'} "
              f"| generic>1e-6: {', '.join(f'{k}={v:.1e}' for k, v in gen.items()) or '-'}")


if __name__ == "__main__":
    main()

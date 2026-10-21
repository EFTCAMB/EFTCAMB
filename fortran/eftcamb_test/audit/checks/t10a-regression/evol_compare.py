"""Track 10a: compare the perturbation time evolution (EFTCAMB get_scale_evolution) of A0 and B1 for a few models.

Runs each model in both trees (IntTolBoost=2 in both, massless or default neutrinos as given) with evol_k, then prints,
per k and per field, max |A0-B1| / max|B1| over a in [a_min, 1] (default a_min = 1e-3, i.e. after recombination).

Usage: python evol_compare.py OUTDIR MODEL [MODEL ...] [--k 0.001,0.01,0.05] [--amin 1e-3] [--nu0]
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("outdir")
    p.add_argument("models", nargs="+")
    p.add_argument("--k", default="0.001,0.01,0.05")
    p.add_argument("--amin", type=float, default=1e-3)
    p.add_argument("--nu0", action="store_true")
    p.add_argument("--fields", default="")
    a = p.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    ks = [float(x) for x in a.k.split(",")]
    ov = {"omnuh2": "0.0"} if a.nu0 else {}
    jobs = [(t, m) for m in a.models for t in ("A0", "B1")]

    def one(j):
        t, m = j
        f = out / f"{t}_{m}.npz"
        if not f.exists():
            r = harness.run_model(t, harness.PARAMS_DIR / f"{m}.ini", f, overrides=ov, evol_k=ks, threads=1,
                                  python_set={"Accuracy.IntTolBoost": 2.0})
            print(t, m, r["status"], flush=True)

    with ThreadPoolExecutor(4) as ex:
        list(ex.map(one, jobs))
    for m in a.models:
        x, y = harness.load(out / f"A0_{m}.npz"), harness.load(out / f"B1_{m}.npz")
        if "evol_a" not in x or "evol_a" not in y:
            print(m, "no evolution output", x["meta"].get("errors"), y["meta"].get("errors"))
            continue
        # the cache fields "k" and "a" overwrite the harness' evol_k / evol_a keys: use the field arrays
        aa = x["evol_a"][0] if x["evol_a"].ndim == 2 else x["evol_a"]
        sel = aa >= a.amin
        fields = sorted(k[5:] for k in x if k.startswith("evol_") and k[5:] not in ("k", "a", "tau", "z") and k in y)
        if a.fields:
            fields = [f for f in fields if f in a.fields.split(",")]
        print(f"== {m}: max|A0-B1|/max|B1| for a >= {a.amin:g}")
        for ik, k in enumerate(ks):
            rows = []
            for f in fields:
                u, v = x[f"evol_{f}"][ik][sel], y[f"evol_{f}"][ik][sel]
                if u.dtype.kind != "f":
                    continue
                s = np.max(np.abs(v)) if np.any(v) else 0.0
                if s == 0 or not np.all(np.isfinite(u)) or not np.all(np.isfinite(v)):
                    continue
                rows.append((float(np.max(np.abs(u - v)) / s), f))
            rows.sort(reverse=True)
            print(f"  k={k:g}: " + "  ".join(f"{f}={e:.1e}" for e, f in rows[:14]))


if __name__ == "__main__":
    main()

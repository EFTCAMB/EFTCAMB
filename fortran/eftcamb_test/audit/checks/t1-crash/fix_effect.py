"""Run a (patched) tree on EFTCAMB test models and compare with the A0 reference runs, to measure the output effect
of a candidate fix (or show bit-identity).

Usage:
    python fix_effect.py TREE OUTDIR [--models a,b | --pattern '5_*.ini'] [--inis x.ini,y.ini] [--jobs 4] [--threads 2]
                         [--ref REFDIR] [--compare-only]
Prints, per model: status, bit_identical, worst Cl rel diff (lensed total TT/EE/BB/TE), tensor BB rel diff, worst P(k),
sigma8 rel diff, and the largest EFT-function difference.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
             "scratchpad/audit")


def tensor_bb(a: dict, b: dict, lmax: int = 1000) -> float:
    ta, tb = a.get("cl_tensor"), b.get("cl_tensor")
    if ta is None or tb is None:
        return float("nan")
    n = min(len(ta), len(tb), lmax + 1)
    x, y = ta[2:n, 2], tb[2:n, 2]
    m = np.abs(y) > 0
    return float(np.max(np.abs(x[m] - y[m]) / np.abs(y[m]))) if m.any() else float("nan")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("tree")
    p.add_argument("outdir")
    p.add_argument("--models", default="")
    p.add_argument("--pattern", default="*.ini")
    p.add_argument("--inis", default="")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--ref", default=str(AUDIT / "runs" / "A0"))
    p.add_argument("--compare-only", action="store_true")
    a = p.parse_args()
    if a.inis:
        inis = [Path(x) for x in a.inis.split(",")]
    else:
        inis = harness.model_inis(a.pattern)
        if a.models:
            want = set(a.models.split(","))
            inis = [i for i in inis if i.stem in want]
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    if not a.compare_only:
        harness.run_batch(a.tree, inis, out, jobs=a.jobs, threads=a.threads, timeout=3000)
    print(f"{'model':34s} {'status':10s} {'bitid':5s} {'worst_cl':>9s} {'tensBB':>9s} {'worst_pk':>9s} {'sigma8':>9s} "
          f"{'eft_max':>9s} eft_field")
    for ini in inis:
        st = out / f"{ini.stem}.json"
        status = json.loads(st.read_text())["status"] if st.exists() else "missing"
        ref = Path(a.ref) / f"{ini.stem}.npz"
        new = out / f"{ini.stem}.npz"
        if status != "ok" or not ref.exists() or not new.exists():
            print(f"{ini.stem:34s} {status:10s} (no comparison: ref={ref.exists()} new={new.exists()})")
            continue
        da, db = harness.load(new), harness.load(ref)
        c = harness.compare(da, db)
        eft = c.get("eft", {}) or {}
        ef, ev = (max(eft.items(), key=lambda kv: kv[1] if np.isfinite(kv[1]) else -1) if eft else ("", 0.0))
        s8 = c.get("sigma8")
        s8v = float(np.max(np.abs(np.atleast_1d(s8)))) if s8 is not None else float("nan")
        print(f"{ini.stem:34s} {status:10s} {str(c['bit_identical'])[:5]:5s} {c['worst_cl']:9.2e} "
              f"{tensor_bb(da, db):9.2e} {c['worst_pk']:9.2e} {s8v:9.2e} {ev:9.2e} {ef}")


if __name__ == "__main__":
    main()

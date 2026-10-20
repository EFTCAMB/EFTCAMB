"""Run a set of test models on a tree with given ini overrides (track t10b-convergence).

python run_set.py LABEL [--tree A0] [--models 'pat1,pat2'] [--set key=value ...] [--jobs 4] [--threads 2] [--skip-done]
Outputs: WORK/runs/LABEL/<model>.npz (+json/log). Models are run longest-first (by default-run elapsed time).
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t10b-convergence"
REF = harness.AUDIT_ROOT / "runs" / "A0"


def ref_elapsed(stem: str) -> float:
    try:
        return json.loads((REF / f"{stem}.json").read_text())["elapsed"]
    except Exception:  # noqa: BLE001
        return 10.0


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("label")
    p.add_argument("--tree", default="A0")
    p.add_argument("--models", default="*")
    p.add_argument("--set", nargs="*", default=[])
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--skip-done", action="store_true")
    p.add_argument("--evol-k", default="")
    p.add_argument("--ini-dir", default="")
    a = p.parse_args()
    pats = a.models.split(",")
    if a.ini_dir:
        inis = sorted(Path(a.ini_dir).glob("*.ini"))
    else:
        inis = harness.model_inis()
    inis = [i for i in inis if any(fnmatch.fnmatch(i.stem, pt) for pt in pats)]
    outdir = WORK / "runs" / a.label
    outdir.mkdir(parents=True, exist_ok=True)
    if a.skip_done:
        inis = [i for i in inis if not (outdir / f"{i.stem}.npz").exists()]
    inis.sort(key=lambda i: -ref_elapsed(i.stem))
    overrides = dict(s.split("=", 1) for s in a.set)
    evk = [float(x) for x in a.evol_k.split(",") if x]
    t0 = time.time()
    res = harness.run_batch(a.tree, inis, outdir, jobs=a.jobs, overrides=overrides, threads=a.threads,
                            evol_k=evk or None, timeout=3000)
    for r in res:
        print(f"{Path(r['model_ini']).stem:40s} {r['status']:20s} {r['elapsed']:8.1f}s", flush=True)
    print(f"total wall {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()

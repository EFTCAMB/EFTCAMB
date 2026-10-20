"""Run several (label, overrides) sets of test models through one shared worker pool (track t10b-convergence).

python run_multi.py --jobs 4 [--models pat1,pat2] [--tree A0] [--threads 1] [--evol-k ...] [--budget SEC]
       --spec LABEL:key=value,key=value --spec LABEL2:...
Outputs WORK/runs/LABEL/<model>.npz. Already finished runs (npz + json) are skipped. Longest tasks first.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t10b-convergence"
REF = harness.AUDIT_ROOT / "runs" / "A0"
COST = {"1": 1.0, "2": 5.0, "3": 15.0}


def ref_elapsed(stem: str) -> float:
    try:
        return json.loads((REF / f"{stem}.json").read_text())["elapsed"]
    except Exception:  # noqa: BLE001
        return 10.0


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--spec", action="append", required=True)
    p.add_argument("--tree", default="A0")
    p.add_argument("--models", default="*")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--budget", type=float, default=1e9, help="do not start new tasks after this many seconds")
    p.add_argument("--evol-k", default="")
    p.add_argument("--ini-dir", default="")
    p.add_argument("--timeout", type=float, default=5400)
    p.add_argument("--env", default="", help="extra environment for the workers, K=V,K2=V2")
    a = p.parse_args()
    pats = a.models.split(",")
    inis = sorted(Path(a.ini_dir).glob("*.ini")) if a.ini_dir else harness.model_inis()
    inis = [i for i in inis if any(fnmatch.fnmatch(i.stem, pt) for pt in pats)]
    evk = [float(x) for x in a.evol_k.split(",") if x] or None
    env = dict(kv.split("=", 1) for kv in a.env.split(",") if kv)
    tasks = []
    for spec in a.spec:
        label, _, kv = spec.partition(":")
        ov = dict(s.split("=", 1) for s in kv.split(",") if s)
        boost = ov.get("accuracy_boost", "1")
        outdir = WORK / "runs" / label
        outdir.mkdir(parents=True, exist_ok=True)
        for ini in inis:
            if (outdir / f"{ini.stem}.npz").exists() and (outdir / f"{ini.stem}.json").exists():
                continue
            cost = ref_elapsed(ini.stem) * COST.get(boost, 5.0)
            tasks.append((cost, label, ini, outdir / f"{ini.stem}.npz", ov))
    tasks.sort(key=lambda t: -t[0])
    print(f"{len(tasks)} tasks", flush=True)
    t0 = time.time()

    def work(t):
        _, label, ini, out, ov = t
        if time.time() - t0 > a.budget:
            return label, ini, {"status": "skipped(budget)", "elapsed": 0}
        return label, ini, harness.run_model(a.tree, ini, out, overrides=ov, threads=a.threads, evol_k=evk,
                                             timeout=a.timeout, extra_env=env)

    with ThreadPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(work, t) for t in tasks]
        for f in as_completed(futs):
            label, ini, r = f.result()
            print(f"{label:10s} {ini.stem:32s} {r['status']:16s} {r.get('elapsed', 0):8.1f}s "
                  f"[{time.time() - t0:7.0f}s]", flush=True)
    print(f"total wall {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()

"""Track t4b-derivs: run deriv_worker.py for test models on a tree (subprocess per model, <=4 concurrent).

    python deriv_run.py --tree A0 --outdir DIR [--pattern '*.ini'] [--jobs 4] [--evol-k 0.05] [--set key=value ...]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

PY = sys.executable


def run_one(tree, ini: Path, outdir: Path, overrides: dict, evol_k: list[float], a_lo: float, a_hi: float,
            n_c: int, hs: list[float], timeout: float = 1500, evol: dict | None = None, name: str | None = None,
            python_set: dict | None = None) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    name = name or ini.stem
    wrapper = harness._wrapper_ini(ini, outdir, name, overrides)
    cfg = {
        "tree": str(harness.tree_path(tree)),
        "ini": str(wrapper),
        "out": str(outdir / f"{name}.npz"),
        "a_lo": a_lo,
        "a_hi": a_hi,
        "n_c": n_c,
        "hs": hs,
        "evol_k": evol_k,
        "python_set": python_set or {},
    }
    cfg.update(evol or {"evol_a_lo": 0.02, "evol_a_hi": 0.9, "evol_n_c": 60, "evol_h": 0.01})
    env = dict(os.environ, OMP_NUM_THREADS="1", PYTHONHASHSEED="0")
    env.pop("PYTHONPATH", None)
    log = outdir / f"{name}.log"
    t0 = time.time()
    try:
        # machine-wide CPU slot (1 thread), nice 10; the worker exits if this process dies
        with harness.slot(1), open(log, "w") as f:
            p = subprocess.run([PY, str(HERE / "deriv_worker.py"), json.dumps(cfg)], cwd=outdir, env=env, stdout=f,
                               stderr=subprocess.STDOUT, timeout=timeout, preexec_fn=harness._nice)
        status = "ok" if p.returncode == 0 else f"rc={p.returncode}"
    except subprocess.TimeoutExpired:
        status = "timeout"
    info = {"model": name, "status": status, "elapsed": round(time.time() - t0, 1), "log": str(log)}
    if status != "ok":
        info["log_tail"] = log.read_text(errors="replace").splitlines()[-15:]
    (outdir / f"{name}.json").write_text(json.dumps(info, indent=1))
    return info


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tree", default="A0")
    p.add_argument("--outdir", required=True)
    p.add_argument("--pattern", default="*.ini")
    p.add_argument("--inis", nargs="*", default=[])
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--evol-k", default="")
    p.add_argument("--a-lo", type=float, default=1e-7)
    p.add_argument("--a-hi", type=float, default=0.92)
    p.add_argument("--n-c", type=int, default=240)
    p.add_argument("--hs", default="0.02,0.01,0.005")
    p.add_argument("--set", nargs="*", default=[])
    a = p.parse_args()
    inis = [Path(x) for x in a.inis] if a.inis else harness.model_inis(a.pattern)
    inis = [i for i in inis if i.stem != "1_EFT_GR"] if not a.inis else inis
    overrides = dict(s.split("=", 1) for s in a.set)
    evk = [float(v) for v in a.evol_k.split(",") if v]
    hs = [float(v) for v in a.hs.split(",")]

    def one(ini):
        return run_one(a.tree, ini, Path(a.outdir), overrides, evk, a.a_lo, a.a_hi, a.n_c, hs)

    with ThreadPoolExecutor(a.jobs) as ex:
        for r in ex.map(one, inis):
            print(f"{r['model']:32s} {r['status']:10s} {r['elapsed']:7.1f}s", flush=True)


if __name__ == "__main__":
    main()

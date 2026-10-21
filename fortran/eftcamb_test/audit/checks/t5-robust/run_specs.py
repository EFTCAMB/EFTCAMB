"""Track 5 (resume 4): run an ad-hoc list of CAMB jobs through harness.run_model, deadline bounded.

Each spec is NAME|TREE|INI|key=val;key=val|ENV=val;ENV=val   (TREE = A0 or a path; INI = model stem or a path;
overrides and extra environment variables optional)
Outputs: $WORK/specs/NAME.npz (skips specs whose .json already exists).

Usage: python run_specs.py SPEC [SPEC ...] [--jobs 3] [--threads 1] [--launch-budget 300] [--hard-budget 540]
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t5-robust" / "specs"


def _kv(text: str) -> dict[str, str]:
    out = {}
    for kv in text.split(";"):
        if kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse(spec: str) -> tuple[str, str, Path, dict[str, str], dict[str, str]]:
    parts = spec.split("|") + ["", ""]
    name, tree, ini = parts[0], parts[1], parts[2]
    inip = Path(ini) if "/" in ini else harness.PARAMS_DIR / f"{ini}.ini"
    return name, tree, inip, _kv(parts[3]), _kv(parts[4])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("specs", nargs="+")
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--launch-budget", type=float, default=300.0)
    p.add_argument("--hard-budget", type=float, default=540.0)
    a = p.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    queue = [parse(s) for s in a.specs]
    queue = [q for q in queue if not (WORK / f"{q[0]}.json").exists()]
    print(f"{len(queue)} pending specs", flush=True)
    lock = threading.Lock()

    def worker() -> None:
        while True:
            with lock:
                if not queue or time.time() - t0 > a.launch_budget:
                    return
                name, tree, ini, ov, env = queue.pop(0)
            remaining = a.hard_budget - (time.time() - t0)
            out = WORK / f"{name}.npz"
            st = harness.run_model(tree, ini, out, overrides=ov, threads=a.threads, timeout=max(30.0, remaining),
                                   extra_env=env)
            note = ""
            if st["status"] == "timeout":
                info = json.loads(out.with_suffix(".json").read_text())
                out.with_suffix(".json").unlink()
                note = f" (deadline after {info['elapsed']:.0f}s, will retry)"
            print(f"{time.strftime('%H:%M:%S')} {name:34s} {st['status']:14s} {st['elapsed']:7.1f}s{note}", flush=True)

    ths = [threading.Thread(target=worker) for _ in range(a.jobs)]
    for t in ths:
        t.start()
    for t in ths:
        t.join()
    print(f"done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

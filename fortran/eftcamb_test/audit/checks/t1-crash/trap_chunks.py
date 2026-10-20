"""Run the trap batch (FP traps + -fcheck tree) over all models in time-bounded chunks, resumable.

Each call runs at most --max models that have no result yet in OUTDIR (jobs <= 4, threads = 1, one CPU slot each),
with reduced accuracy settings (crash hunting, not accuracy) and a per-model timeout, then prints a summary line per
model. Call repeatedly until it reports nothing left.

Usage: python trap_chunks.py TREE OUTDIR [--max 8] [--jobs 4] [--timeout 240] [--only a,b] [--set k=v ...] [--redo]
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import harness  # noqa: E402
import trap_run  # noqa: E402

AUDIT_INIS = HERE.parents[1] / "inis"
REDUCED = {"l_max_scalar": "1200", "l_max_tensor": "500", "k_eta_max_tensor": "1000"}


def all_inis() -> list[Path]:
    return harness.model_inis() + sorted(AUDIT_INIS.glob("00*.ini"))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("tree")
    p.add_argument("outdir")
    p.add_argument("--max", type=int, default=8)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--timeout", type=float, default=240)
    p.add_argument("--only", default="")
    p.add_argument("--set", nargs="*", default=[])
    p.add_argument("--redo", action="store_true")
    p.add_argument("--full", action="store_true", help="do not apply the reduced accuracy overrides")
    a = p.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    inis = all_inis()
    if a.only:
        want = a.only.split(",")
        inis = [i for i in inis if i.stem in want]
    if not a.redo:
        inis = [i for i in inis if not (out / f"{i.stem}.json").exists()]
    todo = inis[: a.max]
    over = {} if a.full else dict(REDUCED)
    over.update(dict(s.split("=", 1) for s in a.set))
    jobs = max(1, min(a.jobs, 4))
    with ThreadPoolExecutor(jobs) as ex:
        res = list(ex.map(lambda i: trap_run.run_model_trap(a.tree, i, out / f"{i.stem}.npz", over, [], a.timeout), todo))
    for r in res:
        print(f"{Path(r['model_ini']).stem:34s} {r['status']:16s} {r['elapsed']:8.1f}s")
    print(f"remaining: {len(inis) - len(todo)}")


if __name__ == "__main__":
    main()

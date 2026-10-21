"""Track 5 (resume 3): deadline-bounded runner for the remaining variation runs.

Runs the planned tasks (skipping those with an existing .json) with at most --jobs concurrent CAMB processes at one
OpenMP thread each, through harness.run_model (machine-wide CPU slots, nice 10).  No new run is started after
--launch-budget seconds; every run gets a timeout that ends before --hard-budget seconds, so one invocation stays
inside a 10 minute tool call.  Runs killed by the chunk deadline have their .json removed so that the next
invocation retries them (with --solo, a single run gets the whole budget).

Usage: python run_chunk.py STAGE[,STAGE...] [--jobs 3] [--launch-budget 330] [--hard-budget 520] [--solo]
       python run_chunk.py --list           (print the pending tasks of every stage)
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402
import run_variations as rv  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t5-robust"
P = harness.PARAMS_DIR
INIS = HERE.parents[1] / "inis"

MINI = rv.MINI
CORE_EXTRA = [m for m in rv.CORE if m not in MINI]
FIVE_MINI = [m for m in MINI if m.startswith("5_")]
FIVE_CORE = [m for m in CORE_EXTRA if m.startswith("5_")]
EXTRA_INIS = [str(INIS / f"{m}.ini") for m in ("007p5-OL-gamma", "007p6-OL-alpha", "007p7-shiftsym",
                                                "008p5-beyond-horndeski")]
IC_MODELS = ["2_PEFT_wDE_wCDM_1", "3_RPH_KT_power_1", "4_Des_mc5e_CPL_1", "4_Des_F_R_LCDM_2", "5_quint_exponential",
             "5_hdsk_TG_2", "2_PEFT_Omega_power_2"]

# stage -> list of (variation, [models])
STAGES: dict[str, list[tuple[str, list[str]]]] = {
    "stt": [("stt_x10", MINI), ("stt_d10", MINI)],
    "pto": [("pto_d10", MINI), ("pto_x10", FIVE_MINI)],
    "grt": [("grt_x10", MINI), ("grt_d10", MINI)],
    "bto": [("bto_x10", MINI), ("bto_d10", MINI)],
    "ab2": [("AB2", [m for m in MINI if not m.startswith("2_PEFT_Omega")])],
    "core-pto": [("pto_d10", CORE_EXTRA), ("pto_x10", FIVE_CORE)],
    "core-grt": [("grt_x10", CORE_EXTRA), ("grt_d10", CORE_EXTRA)],
    "core-bto-stt": [("bto_x10", CORE_EXTRA), ("bto_d10", CORE_EXTRA), ("stt_x10", CORE_EXTRA),
                     ("stt_d10", CORE_EXTRA)],
    "core-ab2": [("AB2", [m for m in CORE_EXTRA if not m.startswith("2_PEFT_Omega")])],
    "toggles": [("dphiF", ["5_hdsk_TG_1"]), ("dphiT", ["5_hdsk_TG_2"]), ("mh", ["5_hdsk_TG_1"]),
                ("ub", ["5_hdsk_TG_1", "5_hdsk_TG_2"])],
    "ic": [("ic2", IC_MODELS), ("ic3", IC_MODELS)],
    "extra": [("pto_d10", EXTRA_INIS), ("pto_x10", EXTRA_INIS), ("grt_x10", EXTRA_INIS),
              ("grt_d10", EXTRA_INIS)],
}


def ini_of(m: str) -> Path:
    return Path(m) if "/" in m else P / f"{m}.ini"


def pending(stages: list[str]) -> list[tuple[str, Path, Path]]:
    tasks = []
    for st in stages:
        for var, models in STAGES[st]:
            for m in models:
                ini = ini_of(m)
                out = WORK / "runs" / var / f"{ini.stem}.npz"
                if out.with_suffix(".json").exists():
                    continue
                if any(t[2] == out for t in tasks):
                    continue
                tasks.append((var, ini, out))
    return tasks


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("stages", nargs="?", default="")
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--launch-budget", type=float, default=330.0)
    p.add_argument("--hard-budget", type=float, default=520.0)
    p.add_argument("--solo", action="store_true", help="run only the first pending task with the full budget")
    p.add_argument("--list", action="store_true")
    p.add_argument("--threads", type=int, default=1)
    a = p.parse_args()
    if a.list:
        for st in STAGES:
            t = pending([st])
            print(f"{st:14s} {len(t):3d} pending: " + " ".join(f"{v}/{i.stem}" for v, i, _ in t)[:400])
        return
    t_start = time.time()
    tasks = pending(a.stages.split(","))
    if a.solo:
        tasks = tasks[:1]
        a.launch_budget = 1.0
    print(f"{len(tasks)} pending tasks", flush=True)
    lock = threading.Lock()
    queue = list(tasks)

    def worker() -> None:
        while True:
            with lock:
                if not queue or time.time() - t_start > a.launch_budget:
                    return
                var, ini, out = queue.pop(0)
            remaining = a.hard_budget - (time.time() - t_start)
            st = harness.run_model("A0", ini, out, overrides=rv.overrides_for(var, ini), threads=a.threads,
                                   timeout=max(30.0, remaining))
            note = ""
            if st["status"] == "timeout":
                # killed by the chunk deadline: forget the run so that the next invocation retries it
                js = out.with_suffix(".json")
                info = json.loads(js.read_text())
                if a.solo or remaining > 500:
                    note = " (genuine timeout kept)"
                else:
                    js.unlink()
                    note = f" (chunk deadline after {info['elapsed']:.0f}s, will retry)"
            print(f"{time.strftime('%H:%M:%S')} {var:8s} {ini.stem:30s} {st['status']:14s} {st['elapsed']:7.1f}s{note}",
                  flush=True)

    threads = [threading.Thread(target=worker) for _ in range(a.jobs)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    left = pending(a.stages.split(","))
    print(f"chunk done in {time.time() - t_start:.0f}s; {len(left)} tasks left in {a.stages}", flush=True)


if __name__ == "__main__":
    main()

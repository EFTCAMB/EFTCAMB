"""Track 5: staged driver (v2, prioritised for an overloaded machine). Restartable: run_variations skips existing outputs.

Before starting, waits for orphaned harness workers (PIDs given with --wait) and writes the .json status file the killed
parent could not write (ok if the npz exists).
Usage: python queue_main2.py [--jobs 3] [--wait PID,PID,...] [--from STAGE]
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

INIS = HERE.parents[1] / "inis"
EXTRA = ",".join(str(INIS / f"{m}.ini") for m in ("007p5-OL-gamma", "007p6-OL-alpha", "007p7-shiftsym",
                                                    "008p5-beyond-horndeski"))
# models whose EFT functions are small at the turn-on (the RGR threshold can move the turn-on) plus two controls
GRT = ("4_Des_F_R_LCDM_2,4_Des_F_R_LCDM_1,4_Des_F_R_wCDM_1,2_PEFT_Omega_power_2,2_PEFT_gamma5_exp_2,"
       "2_PEFT_Omega_exp_1,2_PEFT_wDE_wCDM_1,5_quint_exponential")
STAGES = [
    ("core-AB2", "AB2", ["--list", "core"]),
    ("core-pto", "pto_d10,pto_x10", ["--list", "core"]),
    ("grt", "grt_d10,grt_x10", ["--list", GRT]),
    ("mini-more", "bto_d10,bto_x10,stt_d10,stt_x10", ["--list", "mini"]),
    ("mini-ic", "ic2,ic3", ["--list", "mini"]),
    ("extra", "AB2,pto_d10,pto_x10,ic2", ["--list", EXTRA]),
    ("rep-AB2-pto", "AB2,pto_d10,pto_x10", ["--list", "rep"]),
    ("rep-more", "bto_d10,bto_x10,stt_d10,stt_x10,ic2,ic3", ["--list", "rep"]),
    ("all-AB2-pto", "AB2,pto_d10,pto_x10", []),
]


def wait_orphans(pids):
    while True:
        alive = []
        for p in pids:
            try:
                os.kill(p, 0)
                alive.append(p)
            except ProcessLookupError:
                pass
        if not alive:
            return
        time.sleep(10)


def fix_json(outs):
    for out in outs:
        out = Path(out)
        js = out.with_suffix(".json")
        if js.exists():
            continue
        ok = out.exists()
        info = {"tree": "A0", "model_ini": "", "overrides": {}, "status": "ok" if ok else "error:orphan",
                "elapsed": None, "out": str(out), "log": str(out.with_suffix(".log")), "note": "json written by queue_main2"}
        if not ok:
            info["log_tail"] = out.with_suffix(".log").read_text(errors="replace").splitlines()[-40:]
        js.write_text(json.dumps(info, indent=2))
        print("fixed json", out.name, info["status"], flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--wait", default="")
    p.add_argument("--outs", default="", help="npz outputs of the orphaned workers")
    p.add_argument("--from", dest="start", default="")
    a = p.parse_args()
    if a.wait:
        wait_orphans([int(x) for x in a.wait.split(",")])
        fix_json([x for x in a.outs.split(",") if x])
    started = not a.start
    for name, var, extra in STAGES:
        started = started or name == a.start
        if not started:
            continue
        print(f"### stage {name}: {var}", flush=True)
        subprocess.run([sys.executable, str(HERE / "run_variations.py"), var, "--jobs", str(a.jobs), *extra],
                       check=False)
    print("### ALL STAGES DONE", flush=True)


if __name__ == "__main__":
    main()

"""Track 5: staged driver for the variation runs (representative models first, then all models).

Each stage calls run_variations.py (which skips outputs that already exist), so the driver can be restarted.
Usage: python queue_main.py [--jobs 3] [--from STAGE]
"""

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
INIS = HERE.parents[1] / "inis"
EXTRA = ",".join(str(INIS / f"{m}.ini") for m in ("007p5-OL-gamma", "007p6-OL-alpha", "007p7-shiftsym",
                                                    "008p5-beyond-horndeski"))
CORE_VARS = "AB2,pto_d10,pto_x10,grt_d10,grt_x10"
MORE_VARS = "bto_d10,bto_x10,stt_d10,stt_x10"
STAGES = [
    ("rep-core", CORE_VARS, ["--list", "rep"]),
    ("extra-core", CORE_VARS, ["--list", EXTRA]),
    ("rep-more", MORE_VARS, ["--list", "rep"]),
    ("rep-ic", "ic2,ic3", ["--list", "rep"]),
    ("all-AB2", "AB2", []),
    ("all-core", "pto_d10,pto_x10,grt_d10,grt_x10", []),
    ("extra-more", MORE_VARS + ",ic2,ic3", ["--list", EXTRA]),
    ("all-more", MORE_VARS, []),
    ("all-ic", "ic2,ic3", []),
]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--from", dest="start", default="")
    a = p.parse_args()
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

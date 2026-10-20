"""Compile/run the standalone intent(out) replica (p31_intent_out_*.f90) at several flag sets with the audit gfortran.

Also dumps the optimized GIMPLE of init_bg to show whether the `success = .true.` store survives.
Usage: python p31_intent_out_build.py OUTDIR
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
COMMON = ["-cpp", "-ffree-line-length-none", "-fPIC", "-w"]
VARIANTS = {"O3": ["-O3"], "O2": ["-O2"], "O0": ["-O0"], "O3_FIXED": ["-O3", "-DFIXED"]}

for name, flags in VARIANTS.items():
    d = OUT / name
    d.mkdir(exist_ok=True)
    for src in ("p31_intent_out_mod.f90", "p31_intent_out_main.f90"):
        cmd = ["gfortran", *COMMON, *flags, "-fdump-tree-optimized", "-c", str(HERE / src), "-o", src[:-4] + ".o"]
        subprocess.run(cmd, cwd=d, check=True)
    subprocess.run(["gfortran", "-o", "p31", "p31_intent_out_mod.o", "p31_intent_out_main.o"], cwd=d, check=True)
    run = subprocess.run(["./p31"], cwd=d, capture_output=True, text=True)
    dumps = sorted(d.glob("p31_intent_out_mod.f90.*.optimized"))
    init_bg = ""
    if dumps:
        txt = dumps[0].read_text()
        i = txt.find("init_bg (")
        init_bg = txt[i : txt.find("\n}\n", i) + 2] if i >= 0 else "(init_bg not found)"
    print(f"=== {name} {flags}: {run.stdout.strip()}")
    print("optimized GIMPLE of init_bg:\n" + init_bg)

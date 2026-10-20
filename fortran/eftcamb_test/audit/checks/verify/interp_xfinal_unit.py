"""Independent unit check of the equispaced linear interpolation index near x_final.

Compiles interp_xfinal_unit.f90 against the module files of frozen trees A0 (release) and AD (-fcheck=all), links the
tree's camblib.so, and calls precompute / value / first_derivative at x = log(1 - eps/2) on the RPH grid.
"""

import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
OUT = AUDIT / "work/verify/interp_xfinal"
OUT.mkdir(parents=True, exist_ok=True)

for tree, flags in (("A0", ["-O3"]), ("AD", ["-g", "-O0", "-fcheck=all", "-fbacktrace"])):
    moddir = AUDIT / f"build/{tree}/camb/Releaselib"
    lib = AUDIT / f"trees/{tree}/camb"
    exe = OUT / f"unit_{tree}"
    cmd = ["gfortran", *flags, "-cpp", "-ffree-line-length-none", "-fopenmp", f"-I{moddir}", str(HERE / "interp_xfinal_unit.f90"),
           "-o", str(exe), str(lib / "camblib.so"), f"-Wl,-rpath,{lib}"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=OUT)
    if r.returncode:
        print(tree, "compile failed:", r.stderr[-2000:])
        continue
    for mode in ("none", "precompute", "value", "deriv"):
        r = subprocess.run([str(exe), mode], capture_output=True, text=True, cwd=OUT, env=dict(os.environ, OMP_NUM_THREADS="1"))
        txt = (r.stdout + r.stderr).strip().splitlines()
        if mode == "none":
            print(f"== {tree}: grid/arith")
            print("\n".join(txt))
        else:
            keep = [ln for ln in txt if ("returned" in ln or "runtime error" in ln or "At line" in ln)]
            print(f"== {tree} mode={mode} rc={r.returncode}: " + " | ".join(keep))
sys.exit(0)

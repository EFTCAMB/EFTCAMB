"""Track t4a-param: -fcheck=all / FPE-trap coverage of the cubic and quintic spline parametrizations.

The frozen AD build cannot run them: gfortran 16.1 at -O0 -fcheck=all miscompiles Spline1DValue / SplineQ1DValue
(the callee-saved register x19 is dereferenced as `self` without ever being loaded -> SIGSEGV at the first value()
call).  Work-around: recompile 04p016, 04p019 and the 1D allocator with the AD flags but -O1 into a private module
directory, link them into the driver ahead of the AD camblib.so (the driver's allocate() then uses these objects),
and run the spline cases written by run_param_checks.py --build AD.

usage: python ad_spline_check.py [OPT [FCHECK]]   (default O1 all)
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
W = AUDIT / "work" / "t4a-param" / "adspl"
B = AUDIT / "build" / "AD"
SRC = AUDIT / "trees" / "AD" / "fortran" / "eftcamb"
CASES = AUDIT / "work" / "t4a-param" / "AD" / "cases"
HERE = Path(__file__).resolve().parent
LIB = "/Users/marco/Software/miniforge3/lib"
FLAGS = ["-g", "-fbacktrace", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-finit-real=snan",
         "-finit-integer=-999999", "-finit-logical=true", "-cpp", "-ffree-line-length-none",
         "-fallow-argument-mismatch", "-fopenmp"]


def sh(cmd, **kw):
    p = subprocess.run(["nice", "-n", "10", *map(str, cmd)], capture_output=True, text=True, **kw)
    if p.returncode:
        sys.exit(f"command failed: {' '.join(map(str, cmd))}\n{p.stderr[-2000:]}")
    return p


OPT = sys.argv[1] if len(sys.argv) > 1 else "O1"
FCHECK = sys.argv[2] if len(sys.argv) > 2 else "all"


def build() -> Path:
    od = W / f"{OPT}_{FCHECK.replace(',', '_')}"
    od.mkdir(parents=True, exist_ok=True)
    inc = ["-J" + str(od), "-I" + str(od), f"-I{B}/camb/Releaselib", f"-I{B}/forutils/Releaselib"]
    objs = []
    for rel in ["04f_parametrizations_1D/04p016_spline_parametrizations_1D.f90",
                "04f_parametrizations_1D/04p019_spline5_parametrizations_1D.f90",
                "04p1_parametrizations_1D_allocator.f90"]:
        o = od / (Path(rel).stem + ".o")
        sh(["gfortran", "-c", "-" + OPT, *[f for f in FLAGS if f != "-fcheck=all"], "-fcheck=" + FCHECK, *inc,
            SRC / rel, "-o", o])
        objs.append(o)
    exe = od / "driver_ADspl"
    sh(["gfortran", "-O0", *FLAGS, *inc, HERE / "param_driver.f90", *objs, B / "camb/Releaselib/camblib.so",
        f"-Wl,-rpath,{LIB}", "-o", exe])
    return exe


def main():
    exe = build()
    nm = subprocess.run(["nm", str(exe)], capture_output=True, text=True).stdout
    print("spline value() defined in the executable:", "spline1dvalue" in nm and "splineq1dvalue" in nm)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    for case in sorted(CASES.glob("spline*.case")):
        p = subprocess.run(["nice", "-n", "10", str(exe), str(case)], capture_output=True, text=True, env=env,
                           timeout=300)
        done = "\nDONE" in p.stdout
        msg = [ln.strip() for ln in p.stderr.splitlines() if re.search(r"Fortran runtime error|Error termination|"
                                                                        r"Floating-point exception|At line|"
                                                                        r"Program received signal", ln)]
        bt = [ln.strip() for ln in p.stderr.splitlines() if " at /" in ln or "\tat " in ln][:3]
        tail = [ln.strip() for ln in p.stdout.splitlines() if "error" in ln.lower()][-1:]
        print(f"{case.stem:10s} rc={p.returncode:4d} done={done} {' | '.join(msg[:3])} {' | '.join(tail)}"
              f" {' | '.join(bt)}")


if __name__ == "__main__":
    main()

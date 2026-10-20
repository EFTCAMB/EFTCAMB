"""Verifier (spline5-order3-errorstop): independent unit-level check of spline5_parametrization_1D with
EFTOmega_Spline_Order = 3, bypassing CAMB entirely. Compiles spline5_o3_unit.f90 against the frozen A0 modules and
camblib.so, then (per order 3/5): mode 0 = value/d1/d2 vs finite differences and the candidate cubic d3 = 6*bcoef(5,ip);
mode 3 = call third_derivative; mode 4 = call fourth_derivative.
"""

import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "verify" / "spline5_o3"
TREE = harness.TREES["A0"]
BUILD = harness.AUDIT_ROOT / "build" / "A0"
X6 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    exe = WORK / "spline5_o3_unit"
    lib = TREE / "camb" / "camblib.so"
    cmd = ["gfortran", "-cpp", "-ffree-line-length-none", "-O0", "-fopenmp",
           f"-I{BUILD / 'camb' / 'Releaselib'}", f"-I{BUILD / 'forutils' / 'Releaselib'}",
           str(HERE / "spline5_o3_unit.f90"), str(lib), f"-Wl,-rpath,{lib.parent}", "-o", str(exe)]
    subprocess.run(cmd, check=True, cwd=WORK)
    for order in (5, 3):
        ini = WORK / f"order{order}.ini"
        lines = [f"EFTOmega_Spline_Order = {order}", "EFTOmega_Spline_Pixels = 6"]
        # nodes of Omega(a) = 0.3 a^2 + 0.2 a^3 (non-trivial third derivative)
        lines += [f"EFTOmegax{i + 1} = {x}" for i, x in enumerate(X6)]
        lines += [f"EFTOmegav{i + 1} = {0.3 * x**2 + 0.2 * x**3}" for i, x in enumerate(X6)]
        ini.write_text("\n".join(lines) + "\n")
        for mode in (0, 3, 4):
            env = dict(os.environ, OMP_NUM_THREADS="1")
            with harness.slot(1):
                p = subprocess.run([str(exe), str(ini), str(mode)], capture_output=True, cwd=WORK, env=env, timeout=120)
            out = (p.stdout + p.stderr).decode("latin-1").rstrip().splitlines()
            print(f"--- order={order} mode={mode} rc={p.returncode}")
            for ln in out:
                print("   ", repr(ln) if "\x00" in ln or "\x03" in ln else ln)


if __name__ == "__main__":
    main()

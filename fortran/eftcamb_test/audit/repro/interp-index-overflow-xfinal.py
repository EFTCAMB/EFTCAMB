"""Repro: equispaced linear interpolation reads x(N+1)/y(N+1) when x is within rounding of x_final.

Code: fortran/eftcamb/02_equispaced_interpolation_linear_1D.f90:234 (precompute) and :281/344/407/470/533/596 (value and
derivatives) compute  ind = int((x - x_initial)/grid_width) + 1  with no upper clamp, then read self%x(ind+1) and
self%y(ind+1).  The guard  "x >= x_final"  does not cover x = x_final - (few ulp), for which the quotient rounds up to
N-1 and ind = N.  This is reached with a shipped ini: 3_RPH_KT_power_1 (EFT_mass_stability = T) ->
EFTCAMB_Stability_Check -> compute_additional_derivs (06p2_abstract_EFTCAMB_designer.f90:74) -> dfridr shrinks its step
to ~1e-16 around a = 1 -> RPH compute_H_derivs(temp_a) (007p2_RPH.f90:1887) -> adotoa%precompute(x ~ -1e-16).

Method: run the model on the frozen bounds-checked build (tree AD, gfortran -fcheck=all).  Expected: the stability
check completes.  Observed: Fortran runtime error "Index '6001' of dimension 1 of array 'self%x' above upper bound of
6000".  In release builds the same access silently reads one element past the table (heap garbage) and feeds it into
the mass-stability derivatives.
Also prints, in pure Python, an example x < x_final for which the index formula gives N (same arithmetic as Fortran).

Run: /Users/marco/Software/miniforge3/bin/python interp-index-overflow-xfinal.py
"""

import math
import os
import re
import subprocess
import sys
import tempfile

AUDIT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREE = f"{AUDIT}/trees/AD"
INI = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/3_RPH_KT_power_1.ini"

CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
cp.WantCls = False
cp.WantTransfer = False
camb.get_background(cp)
print('stability check passed', flush=True)
"""


def python_example() -> None:
    # grid as built by EquispacedLinearIntepolateFunction1DInit: x_i = x_initial + (i-1)*grid_width, i = 1..N
    n = 6000
    for x_initial in (math.log(1e-10), math.log(1e-8), math.log(1e-6), -18.0, -20.0):
        x_final = 0.0
        w = (x_final - x_initial) / (n - 1)
        for k in range(1, 200):
            x = x_final - k * 2.0**-60
            ind = int((x - x_initial) / w) + 1
            if x < x_final and ind >= n:
                print(f"x_initial={x_initial:.6f}: x = x_final - {x_final - x:.3e} gives ind = {ind} (N = {n}) "
                      f"-> reads x(N+1), y(N+1)")
                return
    print("no rounding example found for the sampled grids (the Fortran run below is the real test)")


def main() -> None:
    python_example()
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with tempfile.TemporaryDirectory() as d:
        r = subprocess.run([sys.executable, "-c", CHILD, TREE, INI], env=env, capture_output=True, text=True, cwd=d,
                           timeout=1800)
    txt = r.stdout + r.stderr
    m = re.search(r"At line (\d+) of file \S*/(02_equispaced_interpolation_linear_1D\.f90)\s*\n(Fortran runtime error:.*)",
                  txt)
    print("expected: stability check passed")
    if m:
        print(f"observed: {m.group(2)}:{m.group(1)}: {m.group(3)}")
    else:
        print(f"observed: rc={r.returncode}; last lines: {txt.strip().splitlines()[-3:]}")


if __name__ == "__main__":
    main()

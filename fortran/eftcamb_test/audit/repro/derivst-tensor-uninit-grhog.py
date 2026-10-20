"""Repro: equations.f90 derivst (tensor equations) reads the locals grhog_t, grhor_t before they are assigned.

In the EFTCAMB block of derivst (equations.f90:4263-4264)
    EV%eft_cache%gpresdotdotm_t = -(4/3)*EV%eft_cache%Hdot*(grhog_t + grhor_t) + (16/3)*(grhog_t + grhor_t)*adotoa**2 + ...
but grhog_t/grhor_t are only set ~80 lines later (equations.f90:4341-4342). gpresdotdotm_t then feeds
compute_H_derivs (Hdotdot/Hdotdotdot of full-mapping models, 06p1:150-160) and the second order EFT functions.

The repro runs a tensor-only EFT model (2_PEFT_Omega_const_1) on the runtime-checked build (-finit-real=snan) with the
gfortran FP traps switched on (they are inert from Python unless _gfortran_set_fpe is called). Expected: clean run.
Observed: SIGILL (= trapped invalid FP operation on Apple arm64) inside derivst at equations.f90:4264.
With the one-line fix (use EV%eft_cache%grhog_t / grhor_t) the same run completes.

Run: /Users/marco/Software/miniforge3/bin/python derivst-tensor-uninit-grhog.py
"""

import os
import subprocess
import sys

AUDIT = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/t1-crash")
TREES = {"checked build (code under audit)": f"{AUDIT}/chk/AD", "checked build + fix": f"{AUDIT}/chk2/AD"}
INI = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/2_PEFT_Omega_const_1.ini"
LIBGFORTRAN = "/Users/marco/Software/miniforge3/lib/libgfortran.5.dylib"

CHILD = r"""
import ctypes, sys
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
cp.WantScalars = False
cp.WantTensors = True
cp.WantTransfer = False
cp.DoLensing = False
cp.max_l_tensor = 300
cp.max_eta_k_tensor = 600
lib = ctypes.CDLL(sys.argv[3])
lib._gfortran_set_fpe(1 | 4 | 8)   # invalid | zero | overflow, as -ffpe-trap would do in a Fortran main
res = camb.get_results(cp)
lib._gfortran_set_fpe(0)
cl = res.get_tensor_cls(CMB_unit='muK', raw_cl=True)
print('BB(l=2..4) =', cl[2:5, 2])
"""


def main() -> None:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    for label, tree in TREES.items():
        if not os.path.exists(f"{tree}/camb/camblib.so"):
            print(f"{label}: tree {tree} not built, skipped")
            continue
        r = subprocess.run([sys.executable, "-c", CHILD, tree, INI, LIBGFORTRAN], env=env, capture_output=True,
                           text=True, cwd=AUDIT)
        status = "ok" if r.returncode == 0 else f"killed by signal {-r.returncode}" if r.returncode < 0 else \
            f"rc={r.returncode}"
        last = [ln for ln in r.stdout.splitlines() if "BB(" in ln]
        print(f"{label:34s}: expected ok, observed {status}  {last[0] if last else ''}")


if __name__ == "__main__":
    main()

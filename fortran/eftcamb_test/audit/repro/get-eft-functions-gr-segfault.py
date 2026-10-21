"""Repro: EFTCAMB.get_eft_functions() kills the Python interpreter (SIGSEGV) whenever the results object carries no
initialised EFTCAMB parameter cache.

EFTCAMB_GetEFTFunctions (fortran/eftcamb/11_EFTCAMB_python.f90:382-476) never checks EFTflag, allocated(model) or
allocated(par_cache%grhormass): the massive-neutrino loop reads par_cache%grhormass / nu_masses and then calls
EFTCAMB%model procedures (EXC_BAD_ACCESS at address 0x8 in __eftcamb_python_MOD_eftcamb_geteftfunctions).
Two ordinary ways to get there:
  case 1: a GR parameter set (EFTflag = 0, the default of camb.set_params) after camb.get_background;
  case 2: an EFT model (2_PEFT_Omega_const_1.ini) with WantCls = WantTransfer = False after camb.get_results:
          CAMB_GetResults (camb.f90:83-158) then never calls CAMBdata_SetParams, so the EFTCAMB model and its
          parameter cache (results.f90:557-583) are never initialised for this results object.
Expected: a Python exception (or an empty/GR result). Observed: the child process is killed by signal 11.

Run: /Users/marco/Software/miniforge3/bin/python get-eft-functions-gr-segfault.py [TREE]
"""

import os
import subprocess
import sys

TREE = sys.argv[1] if len(sys.argv) > 1 else (
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
    "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
INI = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/2_PEFT_Omega_const_1.ini"

CHILD_GR = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb, numpy as np
pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965)
pars.WantCls = False
res = camb.get_background(pars)
print('EFTflag =', pars.EFTCAMB.EFTflag, '; get_background done', flush=True)
fields, arr = pars.EFTCAMB.get_eft_functions(res, np.array([0.5, 1.0]))
print('returned adotoa =', arr['adotoa'], flush=True)
"""

CHILD_NOOUT = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb, numpy as np
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.WantCls = False
pars.WantTransfer = False
pars.DoLensing = False
res = camb.get_results(pars)
print('EFTflag =', pars.EFTCAMB.EFTflag, '; get_results (no Cls, no transfer) done', flush=True)
fields, arr = pars.EFTCAMB.get_eft_functions(res, np.array([0.5, 1.0]))
print('returned adotoa =', arr['adotoa'], flush=True)
"""


def run(label: str, child: str, *args: str) -> None:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    r = subprocess.run([sys.executable, "-c", child, TREE, *args], env=env, capture_output=True, text=True,
                       timeout=600)
    last = [ln for ln in r.stdout.splitlines() if ln.startswith(("EFTflag", "returned"))]
    obs = "returned normally" if r.returncode == 0 else (
        f"killed by signal {-r.returncode}" if r.returncode < 0 else f"exit code {r.returncode}")
    print(f"{label}: {' | '.join(last)}")
    print(f"   get_eft_functions: expected a Python exception or a result, observed: {obs}")


def main() -> None:
    run("case 1 (GR, get_background)", CHILD_GR)
    run("case 2 (EFT model, get_results with nothing to compute)", CHILD_NOOUT, INI)


if __name__ == "__main__":
    main()

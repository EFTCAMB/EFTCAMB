"""Uninitialised reads in EFT blocks of derivs/derivst, detected with the PS tree (code under audit + inert audit
patches, -O0 -finit-real=snan, no -fcheck) and gfortran FP traps switched on via _gfortran_set_fpe.

Cases (2_PEFT_Omega_const_1, reduced accuracy):
  scal_off  : scalars only, no AUDIT flags        -> expect trap if w_dark_energy_t (unset for EFTflag/=0) is used
  scal_wde  : scalars only, AUDIT_FIX_WDE=1       -> sets w_dark_energy_t=-1 for EFT, expect clean run
  tens_off  : tensors only, no flags              -> expect trap on grhog_t/grhor_t in derivst (t1 finding)
  tens_fix  : tensors only, AUDIT_FIX_TENSOR=1    -> expect clean run
Usage: python snan_wde.py
"""

import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/t8a-equations")
TREE = W / "PS"
INI = Path(__file__).resolve().parents[3] / "parameters" / "2_PEFT_Omega_const_1.ini"
LIBGF = "/Users/marco/Software/miniforge3/lib/libgfortran.5.dylib"

CHILD = r"""
import ctypes, sys
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
mode = sys.argv[4]
cp.WantScalars = mode == "scal"
cp.WantTensors = mode == "tens"
if mode == "scal":
    cp.InitPower.r = 0
cp.WantTransfer = False
cp.DoLensing = False
cp.max_l = 400
cp.max_eta_k = 800
cp.max_l_tensor = 300
cp.max_eta_k_tensor = 600
lib = ctypes.CDLL(sys.argv[3])
lib._gfortran_set_fpe(1 | 4 | 8)
res = camb.get_results(cp)
lib._gfortran_set_fpe(0)
print("finished ok")
"""

res = {}
for tag, mode, env_extra in (("scal_off", "scal", {}), ("scal_wde", "scal", {"AUDIT_FIX_WDE": "1"}),
                             ("tens_off", "tens", {}), ("tens_fix", "tens", {"AUDIT_FIX_TENSOR": "1"})):
    env = dict(os.environ, OMP_NUM_THREADS="1", **env_extra)
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        r = subprocess.run(["nice", "-n", "10", sys.executable, "-c", CHILD, str(TREE), str(INI), LIBGF, mode],
                           env=env, capture_output=True, text=True, cwd=W, timeout=3000)
    st = "ok" if r.returncode == 0 else (f"signal {-r.returncode}" if r.returncode < 0 else f"rc={r.returncode}")
    res[tag] = {"status": st, "tail": (r.stdout + r.stderr)[-300:]}
    print(tag, st)
(W / "snan_wde.json").write_text(json.dumps(res, indent=1))

"""Child for geteft_nocache_lldb.py: reach EFTCAMB_GetEFTFunctions with an uninitialised EFT cache (crash expected).

argv: tree case   case = noout (EFT ini, nothing to compute, get_results) | grini (base_params.ini with EFTflag=0,
get_background, no set_params so the F-001 string getters are not called on -fcheck builds)
"""

import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, sys.argv[1])
import camb  # noqa: E402
import numpy as np  # noqa: E402

P = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/"
case = sys.argv[2]
if case == "noout":
    p = camb.read_ini(P + "2_PEFT_Omega_const_1.ini", no_validate=True)
    p.WantCls = False
    p.WantTransfer = False
    p.DoLensing = False
    res = camb.get_results(p)
else:
    p = camb.read_ini(P + "base_params.ini", no_validate=True)
    p.EFTCAMB.EFTflag = 0
    res = camb.get_background(p)
print("CHILD params EFTflag", p.EFTCAMB.EFTflag, "results EFTflag", res.Params.EFTCAMB.EFTflag, flush=True)
f, r = p.EFTCAMB.get_eft_functions(res, np.array([0.5, 1.0]))
print("CHILD returned adotoa", r["adotoa"].tolist(), flush=True)

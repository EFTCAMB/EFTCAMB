"""Probe: Horava EFT functions returned by get_eft_functions (zero-filled output cache), one a per call."""
import os
import sys

TREE = sys.argv[1]
INI = sys.argv[2]
MODE = sys.argv[3] if len(sys.argv) > 3 else "background"
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, TREE)
sys.path.insert(1, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

with harness.slot(1):
    import numpy as np
    import camb

    p = camb.read_ini(INI, no_validate=True)
    if MODE == "background":
        r = camb.get_background(p)
    else:
        p.WantCls = False
        r = camb.get_results(p)
    for a in [1e-4, 0.1, 1.0]:
        f, arr = p.EFTCAMB.get_eft_functions(r, np.array([a]))
        print(f"a={a:.1e}: adotoa={arr['adotoa'][0]:.5e} Hdot={arr['Hdot'][0]:.5e} EFTc={arr['EFTc'][0]:.5e}"
              f" EFTLambda={arr['EFTLambda'][0]:.5e} EFTcdot={arr['EFTcdot'][0]:.5e}", flush=True)

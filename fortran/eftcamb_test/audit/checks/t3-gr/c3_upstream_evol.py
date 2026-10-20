"""t3-gr: does upstream get_time_evolution work after a transfer-only run (to tell an EFTCAMB bridge crash from an
upstream one)? argv: tree"""
import faulthandler
import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
faulthandler.enable()
sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

sys.path.insert(0, str(harness.tree_path(sys.argv[1])))
import camb  # noqa: E402
import numpy as np  # noqa: E402

with harness.slot(1):
    cp = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.12, As=2e-9, ns=0.965)
    cp.WantCls = False
    cp.WantScalars = len(sys.argv) > 2 and sys.argv[2] == "scalars"
    cp.WantTransfer = True
    cp.set_matter_power(redshifts=[0.0], kmax=10.0)
    r = camb.get_results(cp)
    print("results ok", flush=True)
    ev = r.get_time_evolution([0.17], np.logspace(1, 4, 50), ["delta_cdm", "Weyl"])
    print("get_time_evolution ok", np.shape(ev), flush=True)

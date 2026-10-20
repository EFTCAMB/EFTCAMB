"""t3-gr: the harness worker dies with SIGBUS in a transfer-only (get_scalar_cls=F) run when get_scale_evolution is
requested. Replay the worker steps one by one (one thread, harness CPU slot), with faulthandler.
argv: tree wrapper_ini variant(full|noeftfun|late_a)"""
import faulthandler
import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
faulthandler.enable()
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

tree, ini, variant = sys.argv[1], sys.argv[2], sys.argv[3]  # tree: name or path
sys.path.insert(0, str(harness.tree_path(tree)))
import camb  # noqa: E402
import numpy as np  # noqa: E402

with harness.slot(1):
    cp = camb.read_ini(ini, no_validate=True)
    print("WantCls", cp.WantCls, "WantTransfer", cp.WantTransfer, "EFTflag", cp.EFTCAMB.EFTflag, flush=True)
    r = camb.get_results(cp)
    print("results ok", flush=True)
    if variant != "noeftfun":
        cp.EFTCAMB.get_eft_functions(r, harness.A_GRID)
        print("get_eft_functions ok", flush=True)
    a = np.logspace(-2.5 if variant == "late_a" else -6, 0, 300)
    names, arr = cp.EFTCAMB.get_scale_evolution(r, [0.17], a)
    print("get_scale_evolution ok", flush=True)

import sys
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import os
os.environ["OMP_NUM_THREADS"] = "1"
import harness
with harness.slot(1):
    sys.path.insert(0, str(harness.TREES["A0"]))
    import camb
    camb.set_feedback_level(0)
    cp = camb.read_ini(str(harness.PARAMS_DIR / "5_scg.ini"), no_validate=True)
    cp.EFTCAMB.EFTCAMB_feedback_level = 0
    print("omegas", cp.ombh2, cp.omch2, cp.omnuh2, flush=True)
    for H0 in [float(x) for x in sys.argv[1:]]:
        cp.H0 = H0
        try:
            r = camb.get_background(cp, no_thermo=True)
            print(H0, "ok", r.cosmomc_theta(), flush=True)
        except camb.CAMBError as e:
            print(H0, "CAMBError", str(e).splitlines()[-1], flush=True)

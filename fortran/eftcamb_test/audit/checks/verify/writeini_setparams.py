"""write-ini-drops-eftcamb: same claim through camb.set_params (Python dict path, no ini read)."""

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import harness  # noqa: E402

TREE = str(harness.TREES["A0"])
os.environ["OMP_NUM_THREADS"] = "1"
with harness.slot(1):
    os.nice(10)
    sys.path.insert(0, TREE)
    import numpy as np
    import camb
    from camb.baseconfig import CAMBValueError

    eft = dict(EFTflag=1, PureEFTmodel=1, PureEFTHorndeski=False, EFTwDE=0, PureEFTmodelOmega=1,
               PureEFTmodelGamma1=0, PureEFTmodelGamma2=0, PureEFTmodelGamma3=0, PureEFTmodelGamma4=0,
               PureEFTmodelGamma5=0, PureEFTmodelGamma6=0, EFTOmega0=0.1)
    p = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.12, As=2.1e-9, ns=0.965, lmax=1200, **eft)
    print("camb from", camb.__file__, "| EFTflag =", p.EFTCAMB.EFTflag, "| model =", p.EFTCAMB.model_name())
    with tempfile.TemporaryDirectory() as tmp:
        f = os.path.join(tmp, "w.ini")
        try:
            p.write_ini(f)
            print("validate=True: OK")
        except CAMBValueError as e:
            print("validate=True: CAMBValueError ->", str(e).split(": ", 1)[-1])
        txt = open(f).read()
        print("file exists after exception:", os.path.exists(f), "| EFT keys in file:",
              [ln for ln in txt.splitlines() if ln.split("=")[0].strip() in eft])
        q = camb.read_ini(f)
        print("re-read EFTflag =", q.EFTCAMB.EFTflag)
        c0 = camb.get_results(p).get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"][2:1000, 0]
        c1 = camb.get_results(q).get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"][2:1000, 0]
        print(f"max|TT_reread/TT_orig-1| = {np.max(np.abs(c1 / c0 - 1)):.3e}")

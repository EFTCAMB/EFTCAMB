"""Independent check: build an EFT model purely via set_params (no ini), then re-call set_params(cp=p, H0=...)
and compare TT against (a) a fresh EFT set_params with that H0 and (b) a GR run. Also tests B1 tree."""
import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import harness  # noqa: E402

tree = sys.argv[1]
sys.path.insert(0, tree)
import numpy as np  # noqa: E402

import camb  # noqa: E402

EFT = dict(EFTflag=1, PureEFTmodel=1, PureEFTmodelOmega=1, EFTOmega0=1.0, EFTwDE=0, PureEFTHorndeski=True,
           PureEFTmodelGamma1=0, PureEFTmodelGamma2=0, PureEFTmodelGamma3=0, PureEFTmodelGamma4=0,
           PureEFTmodelGamma5=0, PureEFTmodelGamma6=0, EFTCAMB_stability_time=1e-10)
base = dict(H0=67.0, ombh2=0.022, omch2=0.12, lmax=1200)


def tt(p):
    return camb.get_results(p).get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"][2:1200, 0]


with harness.slot(1):
    p = camb.set_params(**base, **EFT)
    print("fresh EFT:", p.EFTCAMB.EFTflag)
    t_eft = tt(p)
    t_gr = tt(camb.set_params(**base))
    camb.set_params(cp=p, H0=67.0, ombh2=0.022, omch2=0.12)
    print("after set_params(cp=p, H0,ombh2,omch2): EFTflag =", p.EFTCAMB.EFTflag)
    t2 = tt(p)
    print(f"max|TT_after/TT_EFT-1| = {np.max(np.abs(t2 / t_eft - 1)):.3e}")
    print(f"max|TT_after/TT_GR-1|  = {np.max(np.abs(t2 / t_gr - 1)):.3e}")
    print(f"max|TT_EFT/TT_GR-1|    = {np.max(np.abs(t_eft / t_gr - 1)):.3e}")

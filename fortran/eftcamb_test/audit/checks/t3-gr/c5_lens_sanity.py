"""t3-gr (a): EFTCAMB comments out CAMB's lensing AmplitudeError (lensing.f90, 'some theories really destroy the Cls')
and returns early from the lensing routines when Cphil3(10) > lensing_sanity_check_amplitude (1e-7), also for
EFTflag=0. Compare A0 and B0 for GR with As large enough to trip the Fortran check but below the Python
'P(k=0.05) > 2e-8' guard. argv: tree As"""
import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

sys.path.insert(0, str(harness.tree_path(sys.argv[1])))
import camb  # noqa: E402
import numpy as np  # noqa: E402

with harness.slot(1):
    As = float(sys.argv[2])
    cp = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.12, As=As, ns=0.965, lmax=2000, lens_potential_accuracy=0)
    try:
        r = camb.get_results(cp)
        cl = r.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)
        tot, unl = cl["total"], cl["unlensed_total"]
        pp = cl["lens_potential"]
        l = np.arange(pp.shape[0])
        print(f"{sys.argv[1]} As={As:g}: ok; Cphil3(10)~{pp[10, 0] * 2 * np.pi / (10 * 11)**2 * 10**3:.2e}"
              f" lensed/unlensed TT-1 at l=1500: {tot[1500, 0] / unl[1500, 0] - 1:+.3e}; lensed BB(l=1000)={tot[1000, 2]:.3e}")
    except Exception as e:  # noqa: BLE001
        print(f"{sys.argv[1]} As={As:g}: {type(e).__name__}: {str(e)[:150]}")

"""Independent check of horava-hdot-inconsistent: derive Hdot = d(aH)/dtau from CAMB's own background
(results.hubble_parameter, results.conformal_time; not the EFT cache) and compare with eft_Hdot / eft_Hdotdot
stored in the A0 reference npz (harness run on the 401-point eft_a_grid)."""
import os
import sys
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402
import numpy as np  # noqa: E402
from scipy.interpolate import CubicSpline  # noqa: E402

TREE = harness.AUDIT_ROOT / "trees" / "A0"
sys.path.insert(0, str(TREE))
import camb  # noqa: E402

C = 299792.458
print("start", flush=True)
with harness.slot(1):
    for model in ("5_Horava_SolSyst_3", "5_Horava_eta_4", "5_Horava_eta_1", "5_Horava_eta_3"):
        ref = harness.load(harness.AUDIT_ROOT / "runs" / "A0" / f"{model}.npz")
        a_g = ref["eft_a_grid"]
        cp = camb.read_ini(str(REPO / "fortran/eftcamb_test/parameters" / f"{model}.ini"), no_validate=True)
        res = camb.get_background(cp)
        # CAMB background on a dense grid in tau
        a = np.logspace(-3.5, 0, 1500)
        z = 1 / a - 1
        tau = np.asarray(res.conformal_time(z))
        calH = a * np.array(res.hubble_parameter(z)) / C  # conformal Hubble in 1/Mpc
        sH = CubicSpline(tau, calH)
        Hd = sH.derivative(1)
        Hdd = sH.derivative(2)
        sel = (a_g > 2e-3) & (a_g < 0.98)
        tau_g = np.interp(np.log(a_g[sel]), np.log(a), tau)
        Hdot_bg, Hdd_bg = Hd(tau_g), Hdd(tau_g)
        hc = ref["eft_Hdot"][sel]
        hdd = ref["eft_Hdotdot"][sel]
        adot = ref["eft_adotoa"][sel]
        e_ad = np.max(np.abs(adot - sH(tau_g)) / sH(tau_g))
        e1 = np.abs(hc - Hdot_bg) / np.abs(adot**2)
        e2 = np.abs(hdd - Hdd_bg) / np.abs(adot**3)
        print(f"{model}: adotoa cache vs CAMB bg {e_ad:.2e}; |Hdot_cache-dH/dtau|/H^2 max {e1.max():.3e} "
              f"(median {np.median(e1):.3e}); |Hdotdot-d2H/dtau2|/H^3 max {e2.max():.3e}; "
              f"Hdotdotdot max|.| {np.max(np.abs(ref['eft_Hdotdotdot'])):.3e}")

"""Track t4b-derivs: early-time background consistency below/around EFTCAMB_back_turn_on.

For a model ini and a list of EFTCAMB_back_turn_on values, prints (i) CAMB's a H(z)/c (from dtauda) and the
perturbation-side adotoa (get_eft_functions -> compute_adotoa, the value derivs uses) at early a, and (ii) the
conformal time of the EFT perturbation switch-on, tau(a = EFTCAMB_turn_on_time).

    python back_turn_on_probe.py MODEL.ini [bto1 bto2 ...]
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

sys.path.insert(0, str(harness.TREES["A0"]))
import camb  # noqa: E402
import numpy as np  # noqa: E402


def main() -> None:
    ini = Path(sys.argv[1]).resolve()
    btos = sys.argv[2:] or ["1e-8"]
    a = np.array([1e-12, 1e-11, 1e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 1e-6])
    for bto in btos:
        tmp = Path(tempfile.mkdtemp(prefix="t4b_bto_"))
        w = tmp / "w.ini"
        w.write_text(f"DEFAULT({ini})\noutput_root = {tmp}/o\nEFTCAMB_write_background = F\nfeedback_level = 0\n"
                     f"EFTCAMB_feedback_level = 0\nEFTCAMB_back_turn_on = {bto}\n")
        with harness.slot(1):
            cp = camb.read_ini(str(w), no_validate=True)
            res = camb.get_background(cp)
            _, arr = cp.EFTCAMB.get_eft_functions(res, a)
            Hc = np.asarray(res.hubble_parameter(1 / a - 1)) / (camb.constants.c / 1e3) * a
            a_on = float(res.Params.EFTCAMB.EFTCAMB_pert_turn_on)
            tau_on = float(res.conformal_time(1 / a_on - 1))
            tau_rad = a_on / float(res.Params.H0 / (camb.constants.c / 1e3) * np.sqrt(res.get_Omega("photon", 0) +
                                                                                       res.get_Omega("neutrino", 0) +
                                                                                       res.get_Omega("nu", 0)))
        print(f"== {ini.stem}  EFTCAMB_back_turn_on={bto}  pert_turn_on={a_on:.1e}  tau(a_on)={tau_on:.6e} Mpc "
              f"(radiation-era estimate {tau_rad:.6e})")
        for i in range(len(a)):
            print(f"   a={a[i]:.0e}  adotoa(derivs)/aH_CAMB - 1 = {arr['adotoa'][i] / Hc[i] - 1:+.3e}")


if __name__ == "__main__":
    main()

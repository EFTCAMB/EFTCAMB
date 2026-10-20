"""Repro: CAMB v2 default set_for_lmax(lmax) (lens_potential_accuracy=None -> max(4, (lmax-1500)/500)) switches on
non-linear (Halofit) lensing.  With EFTflag != 0 this is applied silently to MG linear spectra, although EFTCAMB treats
non-linearity as not worked out (results.f90 only prints a warning when EFTCAMB_feedback_level > 0) and Halofit gets
EFTCAMB%Effective_w_wa = (-1, 0) for every model except quintessence with EFTCAMB_effective_w0wa=T.
In CAMB 1.x / pre-rebase EFTCAMB the same call (default lens_potential_accuracy=0) kept linear lensing.

Model: 4_Des_F_R_wCDM_1 (designer f(R), w0=-0.9, B0=1).
Run: /Users/marco/Software/miniforge3/bin/python auto-lens-accuracy-enables-nonlinear-for-eft.py
"""

import sys
import tempfile
from pathlib import Path

import numpy as np

import os  # noqa: E402

os.environ["OMP_NUM_THREADS"] = "2"
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
try:  # hold 2 machine-wide CPU slots of the audit harness for the lifetime of this script, when available
    from harness import slot as _slot

    _SLOT = _slot(2)
    _SLOT.__enter__()
except Exception:  # noqa: BLE001
    pass
TREE = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")
PARAMS = Path(__file__).resolve().parents[2] / "parameters"
sys.path.insert(0, TREE)
import camb  # noqa: E402

tmp = Path(tempfile.mkdtemp())
ini = tmp / "fr.ini"
ini.write_text(f"DEFAULT({PARAMS / '4_Des_F_R_wCDM_1.ini'})\nget_tensor_cls = F\nnum_redshiftwindows = 0\n"
               "get_transfer = F\nfeedback_level = 0\nEFTCAMB_write_background = F\n")


def lensed(**kw):
    cp = camb.read_ini(str(ini))
    cp.set_for_lmax(2500, **kw)
    cl = camb.get_results(cp).get_cmb_power_spectra(CMB_unit="muK", raw_cl=True, lmax=2500)["total"]
    return cp, cl


cp_def, cl_def = lensed()  # v2 default call
cp_lin, cl_lin = lensed(nonlinear=False)  # same kmax, linear lensing
cp_old, cl_old = lensed(lens_potential_accuracy=0)  # CAMB 1.x default call
print(f"EFTflag={cp_def.EFTCAMB.EFTflag}, model={cp_def.EFTCAMB.model_name()}, "
      f"w0wa passed to Halofit={cp_def.EFTCAMB.get_effective_w0wa()} (model background w0=-0.9)")
print(f"set_for_lmax(2500) default : NonLinear={cp_def.NonLinear}, max_eta_k={cp_def.max_eta_k:.0f}")
print(f"set_for_lmax(2500, lens_potential_accuracy=0) [1.x default]: NonLinear={cp_old.NonLinear}")
for lab, i in (("TT", 0), ("EE", 1), ("BB", 2)):
    r = np.abs(cl_def[2:, i] / cl_lin[2:, i] - 1)
    print(f"lensed {lab}: max |default/linear - 1| over l<=2500 = {r.max():.3e} (at l={np.argmax(r) + 2})")

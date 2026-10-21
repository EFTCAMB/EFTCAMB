"""Check 18: CAMB v2 non-linear model hooks with EFT on (designer f(R) 4_Des_F_R_LCDM_1, A0).

ExternalNonLinearRatio with a constant ratio sqrt(P_NL/P_L) = 1.1 must give P_NL = 1.21 P_L for the EFT linear P(k)
and lensed Cls consistent with NonLinear_lens; SecondOrderPK (example model, z > 10 only) must run and be finite.
usage: python c18_new_nl_models.py
"""

import json
import os
import sys
import tempfile
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "2"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(harness.tree_path("A0")))
import camb  # noqa: E402
from camb import model  # noqa: E402
from camb.nonlinear import ExternalNonLinearRatio, SecondOrderPK  # noqa: E402

tmp = Path(tempfile.mkdtemp())
ini = tmp / "fr.ini"
MODEL = sys.argv[1] if len(sys.argv) > 1 else '4_Des_F_R_LCDM_1'
ini.write_text(f"DEFAULT({harness.PARAMS_DIR / (MODEL + '.ini')})\nget_tensor_cls = F\nnum_redshiftwindows = 0\n"
               "feedback_level = 0\nEFTCAMB_write_background = F\n")
res = {}
os.nice(10)
with harness.slot(2):
    cp = camb.read_ini(str(ini))
    cp.set_matter_power(redshifts=[20.0, 1.0, 0.0], kmax=5.0, nonlinear=True)
    ext = ExternalNonLinearRatio()
    kh = np.logspace(-4, 2, 50)
    zz = np.array([0.0, 1.0, 5.0, 30.0])
    ext.set_ratio(kh, zz, np.full((len(zz), len(kh)), 1.1))
    cp.NonLinearModel = ext
    cp.NonLinear = model.NonLinear_both
    r = camb.get_results(cp)
    k, z, plin = r.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
    _, _, pnl = r.get_nonlinear_matter_power_spectrum(hubble_units=False, k_hunit=False)
    res["external_ratio_EFTflag"] = int(cp.EFTCAMB.EFTflag)
    res["external_ratio_pnl_over_plin_minus_1.21"] = float(np.max(np.abs(pnl / plin - 1.21)))
    res["external_cl_finite"] = bool(np.all(np.isfinite(r.get_total_cls(2000))))

    cp2 = camb.read_ini(str(ini))
    cp2.set_matter_power(redshifts=[30.0, 20.0, 12.0], kmax=5.0, nonlinear=True)
    cp2.NonLinearModel = SecondOrderPK()
    cp2.NonLinear = model.NonLinear_pk
    cp2.WantCls = False
    r2 = camb.get_results(cp2)
    k2, z2, pl2 = r2.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
    _, _, pn2 = r2.get_nonlinear_matter_power_spectrum(hubble_units=False, k_hunit=False)
    rat = pn2 / pl2
    res["secondorder_finite"] = bool(np.all(np.isfinite(rat)))
    res["secondorder_ratio_range_z"] = {f"{zv:g}": [float(rat[i].min()), float(rat[i].max())] for i, zv in enumerate(z2)}
print(json.dumps(res, indent=1))

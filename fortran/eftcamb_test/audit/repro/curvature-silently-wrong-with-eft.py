"""Repro: omk != 0 is silently accepted with EFTflag != 0, but EFTCAMB's perturbation background ignores curvature.

derivs (equations.f90 ~L2786-2810) takes adotoa from model%compute_adotoa, e.g. EFTCAMBPureEFTstdComputeAdotoa:
adotoa = sqrt((grhom_t + grhov_t)/3) with no State%grhok term (and cothxor = 1/tau, the flat value, at ~L2907).
Model: 2_PEFT_wDE_LCDM_1 = pure-EFT designer with w=-1 and all EFT functions zero, i.e. exactly LCDM (RGR time > 1,
so EFT perturbations are never even switched on).  It must reproduce EFTflag=0 at any omk.
Expected: EFT/GR - 1 ~ 1e-5 (as at omk=0).
Run: /Users/marco/Software/miniforge3/bin/python curvature-silently-wrong-with-eft.py
"""

import sys
import tempfile
from pathlib import Path

import os  # noqa: E402

os.environ["OMP_NUM_THREADS"] = "2"
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
try:  # hold 2 machine-wide CPU slots of the audit harness for the lifetime of this script, when available
    from harness import slot as _slot

    _SLOT = _slot(2)
    _SLOT.__enter__()
except Exception:  # noqa: BLE001
    pass
TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")
PARAMS = Path(__file__).resolve().parents[2] / "parameters"
sys.path.insert(0, TREE)
import camb  # noqa: E402

tmp = Path(tempfile.mkdtemp())


def run(model, omk):
    ini = tmp / f"{model}_{omk}.ini"
    ini.write_text(f"DEFAULT({PARAMS / (model + '.ini')})\nomk = {omk}\nget_tensor_cls = F\nnum_redshiftwindows = 0\n"
                   "feedback_level = 0\nEFTCAMB_write_background = F\ntransfer_kmax = 2\n")
    cp = camb.read_ini(str(ini))
    r = camb.get_results(cp)
    return cp.EFTCAMB.EFTflag, r.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"][:, 0], r.get_sigma8()[-1]


for omk in (0.0, 0.02):
    f0, gr, s_gr = run("1_EFT_GR", omk)
    f1, eft, s_eft = run("2_PEFT_wDE_LCDM_1", omk)
    rel = eft / gr - 1
    print(f"omk={omk:5.2f}  EFTflag {f0} vs {f1}:  TT(EFT)/TT(GR)-1  l=2: {rel[2]:+.3e}  l=10: {rel[10]:+.3e}  "
          f"l=30: {rel[30]:+.3e}  l=1000: {rel[1000]:+.3e} | sigma8 {s_eft:.5f} vs {s_gr:.5f} "
          f"({s_eft / s_gr - 1:+.2e})")

"""Repro: 21cm sources with EFTflag != 0 are refused by the Fortran ini driver but silently computed from Python.

camb.f90 CAMB_RunFromIni:  if (P%EFTCAMB%EFTflag /= 0 .and. P%Do21cm) -> 'EFTCAMB is not compatible with 21 cm
sources so far' + error stop.  The Python paths (camb.read_ini + get_results, or CAMBparams.SourceWindows with a 21cm
window) have no such check, and the 21cm line-of-sight sources in equations.f90 (window_21cm branch, ~L1800-1880)
use GR Einstein-equation substitutions (etak, dgrho, sigma) with no EFT branch.
Expected: an explicit error (as in the ini driver).  Observed: 21cm C_l returned without error or warning.
Run: /Users/marco/Software/miniforge3/bin/python 21cm-not-refused-in-python-with-eft.py
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
from camb.sources import GaussianSourceWindow  # noqa: E402

tmp = Path(tempfile.mkdtemp())
res = {}
for model in ("1_EFT_GR", "2_PEFT_Omega_const_1"):
    ini = tmp / f"{model}.ini"
    # same keys a user would put in an ini: Do21cm = T with a 21cm redshift window
    ini.write_text(f"DEFAULT({PARAMS / (model + '.ini')})\nget_tensor_cls = F\nget_transfer = F\nfeedback_level = 0\n"
                   "EFTCAMB_write_background = F\nDo21cm = T\nnum_redshiftwindows = 1\nredshift(1) = 1.0\n"
                   "redshift_kind(1) = 21cm\nredshift_sigma_Mhz(1) = 5\nline_basic = T\nline_distortions = T\nl_max_scalar = 600\n")
    cp = camb.read_ini(str(ini))
    cp.set_for_lmax(600, lens_potential_accuracy=0)
    cp.SourceTerms.limber_windows = True
    r = camb.get_results(cp)
    cl = r.get_source_cls_dict(raw_cl=True)["W1xW1"]
    res[model] = cl
    print(f"{model:22s} EFTflag={cp.EFTCAMB.EFTflag} Do21cm={cp.Do21cm} window={cp.SourceWindows[0].source_type}: "
          f"21cm C_l(W1xW1) l=10: {cl[10]:.4e}  l=100: {cl[100]:.4e}  l=500: {cl[500]:.4e}  (no error raised)")
r = res["2_PEFT_Omega_const_1"] / res["1_EFT_GR"]
print(f"MG/GR 21cm C_l ratio: l=10 {r[10]:.4f}, l=100 {r[100]:.4f}, l=500 {r[500]:.4f}")

"""Side check for get-evolution-sort-misalign: which fields differ between two IDENTICAL sorted get_evolution calls?
(explains the ~1e-317 'FIX vs current' differences in verify_get_evolution_sort.py as run-to-run garbage, not the fix)."""

import os
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
sys.path.insert(0, f"{ROOT}/trees/A0")
import camb  # noqa: E402

tmp = tempfile.mkdtemp(prefix="rep_", dir=f"{ROOT}/work/verify/getevolsort")
for model in ["5_quint_exponential", "5_hdsk_TG_1"]:
    ini = os.path.join(tmp, f"{model}.ini")
    Path(ini).write_text(f"DEFAULT({PARAMS / (model + '.ini')})\noutput_root = {tmp}/{model}\n"
                         "EFTCAMB_write_background = F\nfeedback_level = 0\n")
    cp = camb.read_ini(ini, no_validate=True)
    res = camb.get_background(cp)
    eta = np.sort(res.conformal_time(np.array([0.0, 0.3, 1.0, 2.5, 6.0, 20.0, 100.0])))
    k = np.array([0.2, 0.01, 0.05])
    _, e1 = cp.EFTCAMB.get_evolution(res, k, eta)
    _, e2 = cp.EFTCAMB.get_evolution(res, k, eta)
    diff = {}
    for f in e1.dtype.names:
        if np.issubdtype(e1.dtype[f], np.floating):
            d = np.abs(e1[f] - e2[f])
            if np.any(d > 0) or np.any(np.isnan(e1[f]) != np.isnan(e2[f])):
                diff[f] = (float(np.nanmax(d)), float(np.nanmax(np.abs(e1[f]))))
    print(model, "fields differing between two identical sorted calls:", diff)

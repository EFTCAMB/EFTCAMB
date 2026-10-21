"""Repro: the quasi-static Sigma output (cache field QS_sigma) is a verbatim copy of the QS_mu formula.

Code: fortran/eftcamb/06_abstract_EFTCAMB_model.f90  EFTCAMBModelComputeQSSigma (line 972) assigns to
eft_cache%QS_sigma exactly the same right-hand side as EFTCAMBModelComputeQSMu (line 920) assigns to QS_mu
(token-for-token identical), so QS_sigma == QS_mu at every (k, a).  The general QS expressions computed by
EFTCAMBModelComputeQSMuSigma give Sigma = mu*(1+gamma)/2 != mu whenever the gravitational slip gamma != 1.

Method: tree A0, get_scale_evolution (Python API) for a designer f(R) model (4_Des_F_R_LCDM_1, gamma != 1);
compare QS_sigma with QS_mu, and QS_sigma_gen with QS_mu_gen.
Usage: /Users/marco/Software/miniforge3/bin/python qs-sigma-copy-of-qs-mu.py [TREE]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0"
)
INI = REPO / "fortran/eftcamb_test/parameters/4_Des_F_R_LCDM_1.ini"
CHILD = r"""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_results(pars)
k = [0.01, 0.1, 1.0]
a = np.array([0.3, 0.5, 0.8, 1.0])
fields, arr = pars.EFTCAMB.get_scale_evolution(res, k, a)
arr = np.asarray(arr)
out = {f: np.asarray(arr[f]).tolist() for f in ("QS_mu", "QS_sigma", "QS_mu_gen", "QS_sigma_gen")}
print("JSON" + json.dumps(out))
"""


def main():
    import numpy as np

    tmp = Path(tempfile.mkdtemp(prefix="qssig_"))
    w = tmp / "w.ini"
    w.write_text(f"DEFAULT({INI})\nfeedback_level = 0\nEFTCAMB_write_background = F\nget_tensor_cls = F\n"
                 f"do_lensing = F\noutput_root = {tmp / 'out'}\n")
    r = subprocess.run([sys.executable, "-c", CHILD, str(TREE), str(w)], capture_output=True, text=True, cwd=tmp,
                       env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900)
    js = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
    if not js:
        print("child failed", r.returncode, r.stdout[-1500:], r.stderr[-1500:])
        sys.exit(1)
    d = {k: np.asarray(v) for k, v in json.loads(js[-1][4:]).items()}
    print(f"tree: {TREE}; model {INI.name}; k = 0.01, 0.1, 1 /Mpc; a = 0.3, 0.5, 0.8, 1.0")
    print("QS_mu    :", np.array2string(d["QS_mu"].ravel(), precision=6))
    print("QS_sigma :", np.array2string(d["QS_sigma"].ravel(), precision=6))
    print("max |QS_sigma - QS_mu| =", float(np.max(np.abs(d["QS_sigma"] - d["QS_mu"]))))
    print("QS_mu_gen    :", np.array2string(d["QS_mu_gen"].ravel(), precision=6))
    print("QS_sigma_gen :", np.array2string(d["QS_sigma_gen"].ravel(), precision=6))
    print("max |QS_sigma_gen - QS_mu_gen| =", float(np.max(np.abs(d["QS_sigma_gen"] - d["QS_mu_gen"]))))
    print("expected: QS_sigma differs from QS_mu where the general formulas differ (QS_sigma ~ QS_sigma_gen)")


if __name__ == "__main__":
    main()

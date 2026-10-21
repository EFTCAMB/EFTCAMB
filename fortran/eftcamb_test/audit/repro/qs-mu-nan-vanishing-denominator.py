"""Repro: the quasi-static outputs QS_mu, QS_sigma (and QS_mu_gen, QS_sigma_gen) are NaN for shipped EFT test models.

EFTCAMBModelComputeQSMu / QSSigma (06_abstract_EFTCAMB_model.f90:908-1010) evaluate a closed-form ratio whose
denominator is a sum of terms each proportional to EFTOmegaP, EFTGamma2V, EFTGamma5V or EFTc; for a model where all of
these vanish (pure EFT with only Gamma1, Gamma3/Gamma4, ... on an LCDM background, i.e. 2_PEFT_gamma1_const_1,
2_PEFT_gamma4_const_1) the numerator vanishes too and the code returns 0/0 = NaN, without a guard.  The same holds for
QS_mu_gen / QS_sigma_gen (EFTCAMBModelComputeQSMuSigma, 06_abstract_EFTCAMB_model.f90:1012-1062) when only Gamma1 is on.
These fields are returned to the user by EFTCAMB.get_scale_evolution / get_evolution (and dumped in debug mode, where
EV%eft_cache%is_nan then stops the run).
Expected: finite values; for these models the quasi-static limit is GR-like (mu = 1, as the code's own 'mu' field
shows; Gamma1 multiplies only pi-dot terms that drop out of the quasi-static limit).

Usage: /Users/marco/Software/miniforge3/bin/python qs-mu-nan-vanishing-denominator.py [TREE]
"""

from __future__ import annotations

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
MODELS = ["2_PEFT_gamma1_const_1", "2_PEFT_gamma4_const_1"]
CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb, numpy as np
p = camb.read_ini(sys.argv[2], no_validate=True)
r = camb.get_results(p)
q = np.array([0.01, 0.1]); a = np.array([0.05, 0.5, 1.0])
names, arr = p.EFTCAMB.get_scale_evolution(r, q, a)
for k in ["QS_mu", "QS_sigma", "QS_mu_gen", "QS_sigma_gen", "mu", "EFTOmegaP", "EFTc", "EFTGamma1V", "EFTGamma4V"]:
    v = np.asarray(arr[k]).ravel()
    print(f"    {k:13s} nan: {int(np.isnan(v).sum())}/{v.size}  values: {np.array2string(v, precision=6)}")
"""


def main():
    tmp = Path(tempfile.mkdtemp(prefix="qsnan_"))
    print(f"tree: {TREE};  k = 0.01, 0.1 /Mpc;  a = 0.05, 0.5, 1")
    for m in MODELS:
        w = tmp / f"{m}.ini"
        w.write_text(f"DEFAULT({REPO / 'fortran/eftcamb_test/parameters' / (m + '.ini')})\nfeedback_level = 0\n"
                     f"EFTCAMB_feedback_level = 0\nEFTCAMB_write_background = F\noutput_root = {tmp / m}\n")
        r = subprocess.run([sys.executable, "-c", CHILD, str(TREE), str(w)], capture_output=True, text=True, cwd=tmp,
                           env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900)
        print(f"  {m} (rc={r.returncode})")
        print(r.stdout.rstrip() or r.stderr[-1500:])
    print("expected: no NaN in QS_mu / QS_sigma / QS_mu_gen / QS_sigma_gen (GR-like value 1 for these models)")


if __name__ == "__main__":
    main()

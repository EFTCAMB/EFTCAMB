"""Check the secondary claim of qs-mu-nan-vanishing-denominator: 5_ADE_1 mu/gamma/sigma_eft non-finite at 1 of 16
points (k = 1e-3, 0.01, 0.1, 1; a = 0.02, 0.1, 0.5, 1), reproducing t2d-static's sampling (get_background or
get_results) and printing the offending denominators dgrho + 3 dgq H/k (+2 dgpi).

Usage: python verify_qsnan_ade.py TREE [background|results]
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

PAR = Path(__file__).resolve().parents[3] / "parameters"
CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb, numpy as np
p = camb.read_ini(sys.argv[2], no_validate=True)
r = camb.get_background(p) if sys.argv[3] == "background" else camb.get_results(p)
q = np.array([1e-3, 0.01, 0.1, 1.0]); a = np.array([0.02, 0.1, 0.5, 1.0])
names, arr = p.EFTCAMB.get_scale_evolution(r, q, a)
np.set_printoptions(precision=6, linewidth=150)
for f in ["mu", "gamma", "sigma_eft", "Psi", "Phi", "dgrho", "dgq", "dgrho_Newt", "QS_mu", "QS_mu_gen"]:
    print(f"{f:11s}", np.asarray(arr[f]))
den = arr["dgrho"] + 3 * arr["dgq"] * arr["adotoa"] / arr["k"]
print("dgrho+3dgqH/k", den)
"""


def main():
    tree = harness.tree_path(sys.argv[1])
    mode = sys.argv[2] if len(sys.argv) > 2 else "background"
    tmp = Path(tempfile.mkdtemp(prefix="qsade_", dir="/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-"
                                "EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/qsnan"))
    w = tmp / "ade.ini"
    w.write_text(f"DEFAULT({PAR / '5_ADE_1.ini'})\nfeedback_level = 0\nEFTCAMB_feedback_level = 0\n"
                 f"EFTCAMB_write_background = F\noutput_root = {tmp / 'ade'}\n")
    import subprocess

    with harness.slot(1):
        r = subprocess.run([sys.executable, "-c", CHILD, str(tree), str(w), mode], capture_output=True, text=True,
                           cwd=tmp, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900)
    print(f"tree={tree} mode={mode} rc={r.returncode}")
    print(r.stdout[-6000:] or r.stderr[-2000:])


if __name__ == "__main__":
    main()

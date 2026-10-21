"""FD check of Hdotdot/Hdotdotdot from get_eft_functions (python bridge path) for a given tree (verify track)."""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
import tempfile
from pathlib import Path
import numpy as np
TREE = Path(sys.argv[1])
sys.path.insert(0, str(TREE)); sys.path.insert(1, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
import camb  # noqa: E402
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
W8 = np.array([1 / 280, -4 / 105, 1 / 5, -4 / 5, 0.0, 4 / 5, -1 / 5, 4 / 105, -1 / 280]); OFFS = np.arange(-4, 5)
with harness.slot(1):
    for model in ("5_Horava_eta_4",):
        tmp = Path(tempfile.mkdtemp(prefix="vf_hddd_"))
        (tmp / "w.ini").write_text(f"DEFAULT({PARAMS / model}.ini)\noutput_root = {tmp}/out\nEFTCAMB_write_background = F\n"
                                   "EFTCAMB_feedback_level = 0\nfeedback_level = 0\n")
        cp = camb.read_ini(str(tmp / "w.ini"), no_validate=True)
        res = camb.get_background(cp)
        a_c = np.array([0.02, 0.1, 0.3, 0.5, 0.9]); h = 0.01
        a = (a_c[:, None] * np.exp(h * OFFS[None, :])).ravel()
        _, arr = cp.EFTCAMB.get_eft_functions(res, a)
        g = {f: arr[f].reshape(len(a_c), 9) for f in arr.dtype.names}
        H = g["adotoa"][:, 4]
        for f, d in (("Hdot", "Hdotdot"), ("Hdotdot", "Hdotdotdot")):
            fd = H * (g[f] @ W8) / h
            print(model, TREE.name, f"{d}: stored", g[d][:, 4], "FD", fd, "rel", np.abs(g[d][:, 4] - fd) / np.abs(fd))

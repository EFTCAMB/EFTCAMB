"""Does the Pade d4 bug reach Cls?  Horndeski (model 7) with Horndeski_designer_V = T, where G4pppp enters
horndeski_designereq_coefficients_p -> phi''' and ddV (008p0_Horndeski.f90:2961-2976), and c3(phi) is a Pade
function (flag 6).  A0 vs one-line fix (work/verify/pade_fix)."""
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify")
TREES = {"A0": str(harness.tree_path("A0")), "FIX": str(W / "pade_fix")}
INI = harness.PARAMS_DIR / "5_hdsk_TG_1.ini"
OUT = W / "pade_hdsk"
D1 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
N1 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.02
OV = {"Horndeski_designer_V": "T", "Horndeski_evolve_hubble": "F", "Horndeski_freefunc0_model": "0",
      "Horndeski_wDE_model": "1", "HorndeskiwDE0": "-0.95", "Horndeski_freefunc3_model": "6",
      "Horndeskic3_Pade_order_N": "1", "Horndeskic3_Pade_order_D": "1", "Horndeskic3a0": "0.0",
      "Horndeskic3N0": "0.0", "Horndeskic3N1": repr(N1), "Horndeskic3D1": repr(D1), "EFTCAMB_feedback_level": "2"}
tag = f"D{D1}_N{N1}"


def one(t):
    return t, harness.run_model(TREES[t], INI, OUT / f"{t}_{tag}.npz", overrides=OV, threads=1)


with ThreadPoolExecutor(2) as ex:
    st = dict(ex.map(one, TREES))
for t, s in st.items():
    print(t, s["status"], s.get("log_tail", [])[-6:] if s["status"] != "ok" else "")
if all(s["status"] == "ok" for s in st.values()):
    a, b = harness.load(OUT / f"A0_{tag}.npz"), harness.load(OUT / f"FIX_{tag}.npz")
    print("errors", a["meta"]["errors"], b["meta"]["errors"])
    c = harness.compare(a, b)
    print(json.dumps({"worst_cl": c["worst_cl"], "worst_pk": c["worst_pk"], "sigma8": c["sigma8"],
                      "bit_identical": c["bit_identical"], "bg_H": c["bg_H"]}, default=str))
    ag = a["eft_a_grid"]
    for f in ["EFTOmegaV", "EFTGamma1V", "EFTGamma2V", "EFTGamma3V", "EFTc", "EFTLambda"]:
        k = "eft_" + f
        if k in a:
            print(f, float(np.max(np.abs(a[k] - b[k])) / max(np.max(np.abs(b[k])), 1e-300)))

"""Adversarial verifier (t3-extgal-gr-point-background): is the small-S region (where the missing z>4 dark energy in
CAMB's dtauda matters most) reachable with the stability check ON?  Runs A0 with fast settings, stability on,
for a few (S, q) points satisfying the model-specific priors, and reports status and CAMB thetastar."""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

INI = REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"
OUT = harness.AUDIT_ROOT / "work/verify/extgal_gr/stab"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "1500", "get_transfer": "F",
        "num_redshiftwindows": "0", "omnuh2": "0.0"}
PTS = {"S0.1_q8": ("0.1", "8"), "S0.5_q1.5": ("0.5", "1.5"), "S2_q0.35": ("2.0", "0.35")}


def job(item):
    n, (s, q) = item
    st = harness.run_model("A0", INI, OUT / f"{n}.npz", overrides={**FAST, "ExtendedGalileon_S": s,
                           "ExtendedGalileon_q": q}, threads=1, timeout=900)
    return n, st


with ThreadPoolExecutor(3) as ex:
    for n, st in ex.map(job, PTS.items()):
        line = {"status": st["status"], "tail": [x.strip() for x in st.get("log_tail", [])[-3:]]}
        if st["status"] == "ok":
            d = harness.load(OUT / f"{n}.npz")
            meta = d["meta"] if isinstance(d["meta"], dict) else json.loads(str(d["meta"]))
            line["thetastar"] = meta.get("derived", {}).get("thetastar")
            import numpy as np
            line["cl_finite"] = bool(np.all(np.isfinite(d["cl_total"][2:1500, :2])))
        print(n, line)

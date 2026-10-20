"""Verifier: run Pure EFT steplog-Omega and designer f(R) steplog-w on A0 and on a tree with only the 04p018 fix."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness

V = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/steplog")
FIX = V / "FIXSL"
cases = {
    "peft": (harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini",
             {"EFT_mass_stability": "F", "PureEFTmodelOmega": "8", "EFTOmega_v1": "0.0", "EFTOmega_v2": "1.5",
              "EFTOmega_at": "0.5", "EFTOmega_delta": "0.5"}),
    "peft_small": (harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini",
             {"EFT_mass_stability": "F", "PureEFTmodelOmega": "8", "EFTOmega_v1": "0.0", "EFTOmega_v2": "0.05",
              "EFTOmega_at": "0.5", "EFTOmega_delta": "0.5"}),
    "desfr": (harness.PARAMS_DIR / "4_Des_F_R_wCDM_1.ini",
              {"EFTwDE": "8", "EFTw_v1": "-1.0", "EFTw_v2": "-0.8", "EFTw_at": "0.3", "EFTw_delta": "0.5"}),
}
which = sys.argv[1:] or list(cases)
for name in which:
    ini, ov = cases[name]
    res = {}
    for tag, tree in (("A0", "A0"), ("FIX", FIX)):
        st = harness.run_model(tree, ini, V / f"{name}_{tag}.npz", overrides=ov, threads=1)
        log = Path(st["log"]).read_text(errors="replace").splitlines() if st.get("log") else []
        res[tag] = st["status"]
        print(name, tag, st["status"], [t.strip() for t in log if "NaN" in t][:2])
    if res["A0"] == "ok" and res["FIX"] == "ok":
        a, b = harness.load(V / f"{name}_A0.npz"), harness.load(V / f"{name}_FIX.npz")
        c = harness.compare(a, b)
        print(name, "worst_cl", c["worst_cl"], "worst_pk", c["worst_pk"], "sigma8", a["sigma8"], b["sigma8"])
        g = "eft_grhov_t"
        if g in a:
            print("grhov_t last 3 A0", a[g][-3:], "FIX", b[g][-3:], "a", a["eft_a_grid"][-3:])
        print(json.dumps({k: v for k, v in c["eft"].items() if v > 1e-6}, default=float)[:800])

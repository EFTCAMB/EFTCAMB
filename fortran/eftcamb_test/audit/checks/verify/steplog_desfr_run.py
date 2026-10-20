"""Run designer f(R) with steplog w(a) (EFTwDE=8) on A0; inspect grhov_t continuity at a=1."""
import sys
import numpy as np
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness

AUD = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
out = AUD + "/work/verify/steplog/desfr_steplog_A0.npz"
ini = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/4_Des_F_R_wCDM_1.ini"
st = harness.run_model(AUD + "/trees/A0", ini, out, overrides={"EFTwDE": 8, "EFTw_v1": -1.0, "EFTw_v2": -0.8,
                       "EFTw_at": 0.3, "EFTw_delta": 0.5}, threads=1)
print(st["status"], st.get("log_tail", "")[-300:] if st["status"] != "ok" else "")
d = harness.load(out)
a = d["eft_a_grid"]
g = d["eft_grhov_t"]
for i in [-4, -3, -2, -1]:
    print(f"a={a[i]:.6f} grhov_t={g[i]:.8e}  grhov_t/a^2-ish ratio to prev {g[i]/g[i-1]:.6f}")
for k in ["eft_adotoa", "eft_Hdot", "eft_EFTOmegaV", "eft_EFTc", "eft_EFTLambda"]:
    if k in d:
        print(k, d[k][-3:])

"""Designer f(R) with the t4a-param steplog w(-1, -0.8, 0.3, 0.5): A0 vs steplog-fix tree (Cl, H(z=0), grhov_t(a=1))."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify")
ov = {"EFTwDE": "8", "EFTw_v1": "-1.0", "EFTw_v2": "-0.8", "EFTw_at": "0.3", "EFTw_delta": "0.5"}
ini = harness.PARAMS_DIR / "4_Des_F_R_wCDM_1.ini"
out = {}
for tag, tree in [("A0", "A0"), ("FIX", W / "steplog_fix")]:
    st = harness.run_model(tree, ini, W / "steplog" / f"C6_{tag}.npz", overrides=ov, threads=1)
    out[tag] = st["status"]
    print(tag, st["status"])
a, b = harness.load(W / "steplog" / "C6_A0.npz"), harness.load(W / "steplog" / "C6_FIX.npz")
c = harness.compare(a, b)
print("worst_cl", c["worst_cl"], "worst_pk", c["worst_pk"], "sigma8", c["sigma8"])
print("cl total", json.dumps(c["cl"]["total"]))
print("H(z=0) A0 / FIX:", a["bg_H"][0], b["bg_H"][0], " grhov_t(a=1)/grhov_t(a=0.955) A0/FIX:",
      a["eft_grhov_t"][-1] / a["eft_grhov_t"][-2], b["eft_grhov_t"][-1] / b["eft_grhov_t"][-2])
print("age A0/FIX", a["meta"]["derived"]["age"], b["meta"]["derived"]["age"])

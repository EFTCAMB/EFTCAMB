"""(i) unstable model -> clean Python exception with a reason; process usable afterwards (same process, stable model)."""
import sys, json
sys.path.insert(0, ".")
import t6lib, harness
ini = harness.PARAMS_DIR / "2_PEFT_gamma1_const_1.ini"
stable = harness.PARAMS_DIR / "2_PEFT_gamma1_const_1.ini"
r = t6lib.probe("A0", ini, overrides={"EFTGamma10": "-1.0", "feedback_level": "0"}, then_ini=stable)
print(json.dumps({k: r.get(k) for k in ("status", "bg_error", "then_ok", "then_H0", "then_error")}))
r = t6lib.probe("A0", ini, overrides={"EFTGamma10": "-1.0", "feedback_level": "3"})
print([l for l in r["log"] if "instab" in l.lower() or "unstable" in l.lower()][:6])

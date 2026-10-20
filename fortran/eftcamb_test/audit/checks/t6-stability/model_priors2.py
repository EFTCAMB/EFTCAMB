"""Follow-up of model_priors.py: (a) designer f(R) B0 sign with EFT_mass_math_stability = T; (b) Horava lambda < -2/3
(allowed by the model prior) with the reason printed by the code (feedback 3)."""
import json
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, ".")
import t6lib  # noqa: E402
import harness  # noqa: E402

H = harness.PARAMS_DIR / "5_Horava_eta_1.ini"
F = harness.PARAMS_DIR / "4_Des_F_R_LCDM_1.ini"
cases = [("f(R) B0, mass_math on", b0, F, {"EFTB0": repr(b0), "EFT_mass_math_stability": "T"}) for b0 in
         (-1e-3, -1e-6, 1e-6, 1e-3)]
cases += [("Horava lambda", lam, H, {"Horava_lambda": repr(lam)}) for lam in (-0.6667, -0.7, -1.0, -5.0)]


def one(c):
    label, x, ini, ov = c
    r = t6lib.probe("A0", ini, overrides={"feedback_level": "3", **ov}, timeout=600)
    reasons = sorted({ln.strip()[:70] for ln in r["log"] if "instab" in ln.lower() or "nan" in ln.lower()
                      or "criteria" in ln.lower()})
    return {"case": label, "x": x, "status": r["status"], "code_stable": r.get("bg_error") is None,
            "reasons": reasons[:4]}


with ThreadPoolExecutor(3) as ex:
    for row in ex.map(one, cases):
        print(json.dumps(row))

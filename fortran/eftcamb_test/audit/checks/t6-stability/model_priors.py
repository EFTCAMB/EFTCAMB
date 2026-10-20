"""Model-specific stability priors vs their analytic boundaries (code verdict through camb.get_background):
Horava (ghost: lambda > 0 or lambda < -2/3; gradient: 0 < eta < 2 xi + 2, Frusciante et al. 2016) and designer f(R)
(B0 > 0 required: B0 < 0 is the tachyonic scalaron; with the shipped generic flags ghost/gradient on, mass off)."""
import json
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, ".")
import t6lib  # noqa: E402
import harness  # noqa: E402

H = harness.PARAMS_DIR / "5_Horava_eta_1.ini"
F = harness.PARAMS_DIR / "4_Des_F_R_LCDM_1.ini"
xi = -1e-5
cases = []
for lam, expect in ((-0.6667, True), (-0.6666, False), (-1e-4, False), (1e-5, True)):
    cases.append(("Horava lambda", lam, H, {"Horava_lambda": repr(lam)}, expect))
for eta, expect in ((-1e-5, False), (1e-5, True), (2 * xi + 2 - 1e-5, True), (2 * xi + 2 + 1e-5, False)):
    cases.append(("Horava eta", eta, H, {"Horava_eta": repr(eta)}, expect))
for b0, expect in ((-1e-3, False), (-1e-6, False), (1e-6, True), (1e-3, True)):
    cases.append(("f(R) B0", b0, F, {"EFTB0": repr(b0)}, expect))
for b0, expect in ((-1e-3, False), (1e-3, True)):
    cases.append(("f(R) B0, mass check on", b0, F, {"EFTB0": repr(b0), "EFT_mass_stability": "T"}, expect))


def one(c):
    label, x, ini, ov, expect = c
    r = t6lib.probe("A0", ini, overrides={"feedback_level": "0", **ov}, timeout=600)
    stable = r["bg_error"] is None if r["status"] == "ok" else None
    return {"case": label, "x": x, "code_stable": stable, "expected_stable": expect, "status": r["status"],
            "err": (r.get("bg_error") or "")[-60:], "agree": stable == expect}


with ThreadPoolExecutor(3) as ex:
    for row in ex.map(one, cases):
        print(json.dumps(row))

"""EFTStabilityComputation (09_EFTCAMB_stability.f90:488-541): the 2nd and 3rd massive-neutrino loops do not reset
eft_cache%gpinu_tot before accumulating (they reset gpinudot_tot twice instead), so compute_H_derivs sees 2x and then
3x the massive-neutrino pressure.  Compare the stability cache (EFTTestStability) with the python bridge
(get_eft_functions, harness reference run) for Pure EFT designer models whose Hdotdot uses gpinu_tot."""
import sys

import numpy as np

sys.path.insert(0, ".")
import t6lib  # noqa: E402
import harness  # noqa: E402

REF = harness.AUDIT_ROOT / "runs" / "A0"
idx = [0, 100, 200, 250, 300, 350, 400]
for model in sys.argv[1:] or ["2_PEFT_wDE_wCDM_1", "2_PEFT_gamma1_const_1"]:
    ref = harness.load(REF / f"{model}.npz")
    a = ref["eft_a_grid"][idx]
    r = t6lib.probe("A0", harness.PARAMS_DIR / f"{model}.ini",
                    overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0"}, a_grid=a, feedback=0)
    print(model, r["status"])
    for j, i in enumerate(idx):
        g_ref = ref["eft_gpinu_tot"][i]
        print(f"  a={a[j]:.2e} gpinu_tot stab={r['gpinu_tot'][j]: .4e} bridge={g_ref: .4e} ratio={r['gpinu_tot'][j] / g_ref if g_ref else float('nan'):.4f}"
              f" | Hdotdot stab={r['Hdotdot'][j]: .5e} bridge={ref['eft_Hdotdot'][i]: .5e} rel={(r['Hdotdot'][j] / ref['eft_Hdotdot'][i] - 1):.2e}"
              f" | Hdotdotdot rel={(r['Hdotdotdot'][j] / ref['eft_Hdotdotdot'][i] - 1):.2e}")

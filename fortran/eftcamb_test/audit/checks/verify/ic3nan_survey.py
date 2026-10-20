"""Survey (no new CAMB runs): for each shipped test model, Omega and Gamma3 at the EFT_IC_type=2/3 evaluation point
a_eval = max(EFTCAMB_turn_on_time, EFTCAMB_back_turn_on), read from the A0 reference npz (get_eft_functions grid), and the
case-3 denominators D1 = 3G3 + 5 Om + 3 G3 Om + 3 Om^2 (equations.f90 and 09 second term) and
D2 = 3 G3 (1 + Om) + Om (5 + 3 Om) (09 first term; identical polynomial to D1).
"""

import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

RUNS = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/runs/A0"
rows = []
for ini in harness.model_inis():
    name = os.path.basename(str(ini))[:-4]
    npz = os.path.join(RUNS, name + ".npz")
    if not os.path.exists(npz):
        rows.append((name, "missing"))
        continue
    p = harness.parse_ini(ini)
    ton = float(p.get("EFTCAMB_turn_on_time", 1e-2))
    bton = float(p.get("EFTCAMB_back_turn_on", 1e-8))
    a_eval = max(ton, bton)
    d = harness.load(npz)
    if "eft_a_grid" not in d:
        rows.append((name, "no-eft-grid (EFTflag=%s)" % p.get("EFTflag")))
        continue
    a = d["eft_a_grid"]
    i = int(np.argmin(abs(np.log(a) - np.log(a_eval))))
    om = float(d["eft_EFTOmegaV"][i])
    g3 = float(d["eft_EFTGamma3V"][i])
    g2 = float(d["eft_EFTGamma2V"][i])
    D1 = 3 * g3 + 5 * om + 3 * g3 * om + 3 * om**2
    rows.append((name, a_eval, a[i], om, g3, g2, D1))

zero = [r for r in rows if len(r) > 2 and r[3] == 0.0 and r[4] == 0.0]
for r in rows:
    if len(r) == 2:
        print(f"{r[0]:32s} {r[1]}")
    else:
        flag = "ZERO-DEN" if r[6] == 0.0 else ("small" if abs(r[6]) < 1e-8 else "")
        print(f"{r[0]:32s} a_eval={r[1]:.3g} (grid {r[2]:.3g}) Omega={r[3]: .3e} Gamma3={r[4]: .3e} Gamma2={r[5]: .3e} "
              f"D={r[6]: .3e} {flag}")
print(f"\n{len(zero)} / {len([r for r in rows if len(r) > 2])} shipped models have Omega = Gamma3 = 0 exactly at a_eval:")
print(", ".join(r[0] for r in zero))
json.dump({"zero_den_models": [r[0] for r in zero], "rows": [list(r) for r in rows]},
          open("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
               "scratchpad/audit/work/verify/ic3nan/survey.json", "w"), indent=1)

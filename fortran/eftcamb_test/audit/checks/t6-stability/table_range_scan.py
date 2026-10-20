"""Search for RPH / Pure EFT parameter points where the verdict changes when the designer background table is
extended below the default model_background_a_ini (1e-8) down past EFTCAMB_stability_time (1e-10)."""
import sys, json, itertools
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0, ".")
import t6lib, harness
from boundary_scans import INIS
ini = INIS / "t6_RPH_KB_const.ini"
pts = list(itertools.product(np.round(np.arange(-0.3, 0.61, 0.05), 3), (-0.02, 0.01, 0.1, 1.0), (0.0, -0.5, 0.5),
                             (0.1, -0.3)))
def one(p):
    B, K, T, M = p
    out = []
    for ainit in (None, "1e-12"):
        ov = {"feedback_level": "0", "RPHkineticity0": repr(K), "RPHbraiding0": repr(float(B)), "RPHtensor0": repr(T),
              "RPHmassP0": repr(M)}
        if ainit:
            ov["model_background_a_ini"] = ainit
        r = t6lib.probe("A0", ini, overrides=ov, name="trs")
        out.append(None if r["status"] != "ok" else r["bg_error"] is None)
    return [float(B), K, T, M] + out
with ThreadPoolExecutor(4) as ex:
    res = list(ex.map(one, pts))
diff = [r for r in res if r[4] != r[5]]
print(json.dumps({"n": len(res), "n_stable_default": sum(1 for r in res if r[4]), "n_differ": len(diff), "differ": diff}))
(t6lib.WORK / "table_range_scan.json").write_text(json.dumps(res))

"""Track t4b-derivs: classify failing derivative relations of deriv_analyze.py (active window).

For each failing (model, stored-field) pair report whether the stored field is identically zero in the
get_eft_functions output (field not filled by the python bridge), constant, or nonzero, and the ratio
stored/predicted at the worst point.

    python deriv_classify.py ANALYSIS.json NPZ_DIR
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from deriv_analyze import RELATIONS, load, predict  # noqa: E402

data = json.load(open(sys.argv[1]))
npz_dir = Path(sys.argv[2])
kinds = {s: (b, k, o) for s, b, k, o in RELATIONS}
by_field = defaultdict(list)
for m in data:
    fails = [s for s, v in m["relations"].items() if v.get("active") and not v["active"]["pass"]]
    if not fails:
        continue
    d = load(npz_dir / f"{m['model']}.npz")
    a = np.exp(d["xc"])
    for s in fails:
        base, kind, order = kinds[s]
        S = d[f"st_{s}"][0, 4]
        pred, _ = predict(d, base, kind)
        ia = m["relations"][s]["active"]["a_max_excess_phys"]
        i = int(np.argmin(np.abs(a - ia)))
        if np.all(S == 0):
            cls = "ZERO(not filled)"
        else:
            cls = f"ratio={S[i] / pred[1][i] if pred[1][i] else np.inf:+.3e}"
        by_field[s].append((m["model"], m["relations"][s]["active"]["max_excess_phys"], ia, cls))
for s, lst in by_field.items():
    print(f"== {s}")
    for model, e, ia, cls in sorted(lst, key=lambda x: -x[1]):
        print(f"   {model:28s} err={e:.2e} @a={ia:.2e}  {cls}")

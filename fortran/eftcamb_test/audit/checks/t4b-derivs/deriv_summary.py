"""Summarise deriv_analyze.py json: per relation, the worst model and the list of failing models.

    python deriv_summary.py ANALYSIS.json [--local]
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict

data = json.load(open(sys.argv[1]))
ACT = "--active" in sys.argv
key = "max_excess" if "--local" in sys.argv else "max_excess_phys"
akey = "a_max_excess" if "--local" in sys.argv else "a_max_excess_phys"
pkey = "pass_local" if "--local" in sys.argv else "pass"
rel = defaultdict(list)
gen = defaultdict(list)
for m in data:
    for k, v in m["relations"].items():
        vv = v["active"] if ACT else v
        if vv is None:
            continue
        rel[k].append((vv[key], vv[akey], m["model"], vv[pkey], v["tol"]))
    for k, v in (m["generic_active"] if ACT else m["generic"]).items():
        gen[k].append((v["max_err"], v["a_at_max"], m["model"], v["pass"], 1e-6))
for title, dd in (("DERIVATIVE RELATIONS", rel), ("GENERIC RELATIONS", gen)):
    print(f"==== {title} ({key if title.startswith('DER') else 'max_err'})")
    for k, lst in dd.items():
        lst.sort(reverse=True)
        fails = [x for x in lst if not x[3]]
        worst = lst[0]
        print(f"{k:20s} n={len(lst):2d} fail={len(fails):2d} tol={worst[4]:.0e} worst={worst[0]:.2e} @a={worst[1]:.2e} "
              f"[{worst[2]}]")
        for x in fails[:40]:
            print(f"      {x[0]:.2e} @a={x[1]:.2e} {x[2]}")

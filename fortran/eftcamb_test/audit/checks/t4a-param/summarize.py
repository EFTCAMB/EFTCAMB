"""Track t4a-param: worst error per method per parametrization family, from the outputs of run_param_checks.py
(finite differences, node/continuity checks, integral vs quad) and exact_checks.py (sympy/mpmath).

usage: python summarize.py [--build A0]
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
METHODS = ["value", "d1", "d2", "d3", "d4"]


def family(name: str) -> str:
    return re.sub(r"_\d+$", "", name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default="A0")
    args = ap.parse_args()
    w = AUDIT / "work" / "t4a-param" / args.build
    res = json.loads((w / "results.json").read_text())
    exact = json.loads((w / "exact_results.json").read_text())
    fam: dict = defaultdict(lambda: {"rel": defaultdict(float), "sc": defaultdict(float), "fd": defaultdict(float),
                                     "int": 0.0, "nonfinite": 0, "node": 0.0, "jump": 0.0, "cases": 0})
    for r in res:
        f = fam[family(r["name"])]
        f["cases"] += 1
        for m, v in (r.get("worst") or {}).items():
            if v[0] is not None:
                f["fd"][m] = max(f["fd"][m], float(v[0]))
        wi = r.get("worst_integral")
        if wi and wi[0] is not None:
            f["int"] = max(f["int"], float(wi[0]))
        for iss in r["issues"]:
            mm = re.search(r"does not reproduce node values: worst rel ([0-9.e+-]+)", iss)
            if mm:
                f["node"] = max(f["node"], float(mm.group(1)))
            mm = re.search(r"discontinuity of value at interior nodes: max rel jump ([0-9.e+-]+)", iss)
            if mm:
                f["jump"] = max(f["jump"], float(mm.group(1)))
        e = exact.get(r["name"]) or {}
        for m in METHODS:
            if m in e and isinstance(e[m], dict):
                f["rel"][m] = max(f["rel"][m], float(e[m]["rel"][0]))
                f["sc"][m] = max(f["sc"][m], float(e[m]["scaled"][0]))
                f["nonfinite"] += e[m]["n_nonfinite"]
    out = {}
    for k, f in fam.items():
        row = {"cases": f["cases"], "exact_rel": dict(f["rel"]), "exact_scaled": dict(f["sc"]), "fd": dict(f["fd"]),
               "integral_rel": f["int"], "nonfinite": f["nonfinite"], "node_err": f["node"], "node_jump": f["jump"]}
        out[k] = row
        ex = " ".join(f"{m}={f['sc'].get(m, float('nan')):.1e}/{f['rel'].get(m, float('nan')):.1e}" for m in METHODS
                      if m in f["sc"])
        fd = " ".join(f"{m}={v:.1e}" for m, v in f["fd"].items())
        print(f"{k:14s} n={f['cases']} exact(scaled/rel)[{ex}] fd[{fd}] int={f['int']:.1e} nonfin={f['nonfinite']}"
              f" node={f['node']:.1e} jump={f['jump']:.1e}")
    (w / "summary.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

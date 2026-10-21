"""Track t4b-derivs: derivative relations of the timestep cache at the EFT perturbation turn-on (RGR) time.

The initial conditions of the EFT perturbations (09_EFTCAMB_IC.f90) are set from the cache at a_on =
EFTCAMB_pert_turn_on (after RGR).  For each npz (deriv_worker.py output) the relation errors of
deriv_analyze.check_relation are evaluated on the stencil centres with |ln a - ln a_on| <= W (default 0.15),
using the physically-scaled error (max_excess_phys).

    python deriv_turnon.py DIR [DIR ...] [--w 0.15]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from deriv_analyze import RELATIONS, check_relation, generic_relations, load  # noqa: E402

args = [x for x in sys.argv[1:] if not x.startswith("--")]
w = float(sys.argv[sys.argv.index("--w") + 1]) if "--w" in sys.argv else 0.15
args = [x for x in args if x != str(w)]
seen = set()
worst = {}
for dd in args:
    for f in sorted(Path(dd).glob("*.npz")):
        if f.stem in seen:
            continue
        d = load(f)
        a_on = d["meta"]["pert_turn_on_results"]
        xc = d["xc"]
        if not (xc.min() + 0.05 < np.log(a_on) < xc.max() - 0.05):
            continue
        seen.add(f.stem)
        mask = np.abs(xc - np.log(a_on)) <= w
        fails = []
        for stored, base, kind, order in RELATIONS:
            r = check_relation(d, stored, base, kind, order, mask=mask)
            if r is None:
                continue
            S = d[f"st_{stored}"][0, 4]
            if np.all(S[mask] == 0):
                continue  # not filled by get_eft_functions (reported separately)
            e = r["max_excess_phys"]
            if e > worst.get(stored, (0, ""))[0]:
                worst[stored] = (e, f.stem)
            if not r["pass"]:
                fails.append(f"{stored}={e:.1e}")
        gen = generic_relations(d, mask=mask)
        gfails = [f"{k}={v['max_err']:.1e}" for k, v in gen.items() if not v["pass"]]
        print(f"{f.stem:28s} a_on={a_on:.3e}  FAIL: {' '.join(fails + gfails) or '-'}")
print("\nworst per relation at turn-on:")
for k, (e, m) in sorted(worst.items(), key=lambda x: -x[1][0]):
    print(f"   {k:16s} {e:.2e} [{m}]")

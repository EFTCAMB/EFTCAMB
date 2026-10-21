"""Evaluate theta(H0) at the set_cosmology bracket ends (max(10, 1.1*100*sqrt(omega_m)), 100) and at the ini H0 for
all 73 test models on A0 (fresh read_ini per evaluation; process deaths recorded).  Counts how many test models
cannot be inverted with set_cosmology(cosmomc_theta) at their fiducial parameters.

Usage: python theta_bracket_allmodels.py <outdir> [tree]
"""

import glob
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theta_bracket_indep as tbi  # noqa: E402

outdir = sys.argv[1]
tree = sys.argv[2] if len(sys.argv) > 2 else "A0"
inis = sorted(glob.glob(os.path.join(tbi.PARAMS, "*.ini")))
inis = [i for i in inis if "base_params" not in os.path.basename(i)]


def job(ini):
    keys = tbi.harness.parse_ini(ini)
    om = float(keys["ombh2"]) + float(keys["omch2"]) + float(keys.get("omnuh2", 0))
    lo = max(10.0, 1.1 * 100 * np.sqrt(om))
    H0 = float(keys.get("hubble", 70))
    return os.path.basename(ini)[:-4], tbi.scan(os.path.join(tbi.TREES, tree), ini, [lo, H0, 100.0])


out = {}
with ThreadPoolExecutor(4) as ex:
    for m, r in ex.map(job, inis):
        out[m] = r
json.dump(out, open(os.path.join(outdir, f"allmodels_ends_{tree}.json"), "w"), indent=1)
bad = {}
for m, r in out.items():
    st = [row[1] if row[1] == "ok" else (row[1] + ":" + str(row[2])) for row in r]
    if any(s != "ok" for s in st):
        bad[m] = st
print(f"{len(out)} models; {len(bad)} with a failing evaluation at lo / H0_ini / 100:")
for m, st in bad.items():
    print(f"  {m:40s} lo={st[0]} | H0_ini={st[1]} | 100={st[2]}")

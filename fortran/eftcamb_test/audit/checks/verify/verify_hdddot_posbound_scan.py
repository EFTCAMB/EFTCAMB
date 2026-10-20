"""Verifier: does the Hdotdotdot defect change the stability decision of the positivity bounds / Minkowski-limit check?
Scan the Omega-Lambda gamma model (007p5, uses the base full-map Hdotdotdot, Gamma3 != 0) over (Omega0, Gamma30) with
EFT_positivity_bounds=T or EFT_minkowski_limit=T and print stable/unstable. Run once per tree, compare outputs.
Usage: python verify_hdddot_posbound_scan.py --tree PATH --out file.json"""

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--tree", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, args.tree)
import camb  # noqa: E402

BASE = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/base_params.ini"
tmp = tempfile.mkdtemp(prefix="pbs_", dir=str(Path(args.out).parent))
res = {}
for flag in ("EFT_positivity_bounds", "EFT_minkowski_limit"):
    for om in (0.0, 0.02, 0.05, 0.1):
        for g3 in np.round(np.linspace(-0.3, 0.3, 25), 4):
            for g1, g2 in ((0.0, 0.0), (1.0, -0.1)):
                ini = Path(tmp) / "m.ini"
                ini.write_text(f"""EFTflag = 2
AltParEFTmodel = 2
OLLambdamodel = 0
OLOmegamodel  = {3 if om else 0}
OLGamma1model = {1 if g1 else 0}
OLGamma2model = {3 if g2 else 0}
OLGamma3model = 3
OLOmega0 = {om}
OLOmegaExp = 1.0
OLGamma10 = {g1}
OLGamma20 = {g2}
OLGamma2Exp = 1.0
OLGamma30 = {g3}
OLGamma3Exp = 1.0
EFT_ghost_math_stability = F
EFT_mass_math_stability = F
EFT_ghost_stability = F
EFT_gradient_stability = F
EFT_mass_stability = F
EFT_additional_priors = F
{flag} = T
EFTCAMB_stability_time = 1.e-10
EFTCAMB_stability_threshold = 0.0
EFTCAMB_feedback_level = 0
feedback_level = 0
output_root = {tmp}/m
EFTCAMB_write_background = F
DEFAULT({BASE})
""")
                key = f"{flag}|Om={om}|g1={g1}|g2={g2}|g3={g3}"
                try:
                    cp = camb.read_ini(str(ini), no_validate=True)
                    camb.get_background(cp)
                    res[key] = "stable"
                except Exception as e:  # noqa: BLE001
                    res[key] = "unstable:" + type(e).__name__
Path(args.out).write_text(json.dumps(res, indent=1))
print(args.tree, "stable:", sum(v == "stable" for v in res.values()), "of", len(res))

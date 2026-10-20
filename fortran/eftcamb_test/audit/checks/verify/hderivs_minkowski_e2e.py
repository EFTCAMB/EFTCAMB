"""End-to-end EFT_minkowski_limit verdict from camb.get_background, A0 vs the P_ldd copy (A0 + EFTLambdadotdot supplied
by a central difference before compute_positivity_bounds; see hderivs_build_ldd_patch.py).

For each quintessence test model: verdict with EFT_minkowski_limit = F (control) and = T, EFT_mass_stability = F,
EFTCAMB_stability_time = 1e-3 and the ini default, EFTCAMB_feedback_level = 2 (the stability messages are counted).

Usage: python hderivs_minkowski_e2e.py
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ROOT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
            "scratchpad/audit")
PAR = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")
TREES = {"A0": ROOT / "trees/A0", "P_ldd": ROOT / "work/verify/hderivs/P_ldd"}
MODELS = ["5_quint_exponential", "5_quint_inverse_power", "5_quint_quadratic", "5_quint_axion", "5_quint_cubic",
          "5_Kmouflage_1"]

WORKER = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
ini, mink, stime = sys.argv[2], sys.argv[3] == "T", sys.argv[4]
cp = camb.read_ini(ini, no_validate=True)
cp.EFTCAMB.EFT_minkowski_limit = mink
cp.EFTCAMB.EFT_mass_stability = False
cp.EFTCAMB.EFTCAMB_feedback_level = 2
if stime != "default":
    cp.EFTCAMB.EFTCAMB_stability_time = float(stime)
try:
    camb.get_background(cp)
    v = "stable"
except Exception as e:
    v = "unstable: " + str(e)[:80]
sys.stdout.flush()
print("\nRESULT " + json.dumps({"verdict": v, "camb": camb.__file__}), flush=True)
"""

rows = []
with tempfile.TemporaryDirectory() as d:
    for m in MODELS:
        ini = Path(d) / f"{m}.ini"
        ini.write_text(f"DEFAULT({PAR}/{m}.ini)\nfeedback_level = 0\nEFTCAMB_write_background = F\n"
                       f"output_root = {d}/out\n")
        for stime in ("1e-3", "default"):
            for mink in ("F", "T"):
                for tname, tree in TREES.items():
                    env = dict(os.environ, OMP_NUM_THREADS="1")
                    env.pop("PYTHONPATH", None)
                    with harness.slot(1):
                        p = subprocess.run([sys.executable, "-c", WORKER, str(tree), str(ini), mink, stime],
                                           capture_output=True, text=True, cwd=d, env=env,
                                           preexec_fn=lambda: os.nice(10), timeout=600)
                    res = [json.loads(x[7:]) for x in p.stdout.splitlines() if x.startswith("RESULT ")]
                    nmink = p.stdout.count("Healthy Minkowski limit does not exist")
                    v = res[0]["verdict"] if res else f"crash rc={p.returncode}: {p.stderr[-200:]}"
                    rows.append({"model": m, "stability_time": stime, "minkowski": mink, "tree": tname, "verdict": v,
                                 "n_minkowski_msgs": nmink})
                    print(f"{m:22s} stab_time={stime:7s} minkowski={mink} {tname:6s}: {v[:60]:60s}"
                          f" ('Healthy Minkowski limit does not exist' x{nmink})", flush=True)
out = ROOT / "work/verify/hderivs/minkowski_e2e.json"
out.write_text(json.dumps(rows, indent=1))

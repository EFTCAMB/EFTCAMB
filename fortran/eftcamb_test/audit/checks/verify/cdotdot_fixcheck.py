"""Verifier for stability-cache-cdotdot-uninitialised: end-to-end stability verdict (camb.get_background with the
stability check forced on) on A0 (release), t6 nanini (A0 source, -g -O1 -finit-real=nan -finit-derived) and t6 nanfix
(same flags + the 2-line init fix).  Prints the verdict and the first 'is NaN' message per model/tree.
usage: python cdotdot_fixcheck.py [models...]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

AUDIT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREES = {"A0": f"{AUDIT}/trees/A0", "nanini": f"{AUDIT}/work/t6-stability/nanini",
         "nanfix": f"{AUDIT}/work/t6-stability/nanfix"}
PAR = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"
MODELS = sys.argv[1:] or ["5_ADE_1", "5_scg", "5_Horava_eta_1", "5_quint_axion", "5_Kmimic_1", "4_Des_F_R_wCDM_1",
                          "4_Des_mc5e_wCDM_1", "3_RPH_MK_const_1", "2_PEFT_Omega_exp_1"]

WORKER = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
assert camb.__file__.startswith(sys.argv[1]), camb.__file__
cp = camb.read_ini(sys.argv[2], no_validate=True)
cp.EFTCAMB.EFTCAMB_skip_stability = False
try:
    camb.get_background(cp)
    v = "stable"
except Exception as e:
    v = "UNSTABLE: " + str(e).splitlines()[-1]
print("VERDICT " + json.dumps(v))
"""


def run(tree: str, model: str) -> tuple[str, str | None]:
    with tempfile.TemporaryDirectory() as d:
        ini = Path(d) / "m.ini"
        ini.write_text(f"DEFAULT({PAR}/{model}.ini)\nEFTCAMB_feedback_level = 3\nEFTCAMB_write_background = F\n"
                       f"output_root = {d}/out\n")
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        with harness.slot(1):
            p = subprocess.run(["nice", "-n", "10", sys.executable, "-c", WORKER, tree, str(ini)],
                               capture_output=True, text=True, cwd=d, env=env, timeout=600)
    v = [json.loads(l[8:]) for l in p.stdout.splitlines() if l.startswith("VERDICT ")]
    nan = [l.strip() for l in p.stdout.splitlines() if l.strip().endswith("is NaN")]
    return (v[0] if v else f"NO VERDICT rc={p.returncode} {p.stderr[-300:]}"), (nan[0] if nan else None)


def main():
    out = {}
    for model in MODELS:
        row = {name: run(t, model) for name, t in TREES.items()}
        out[model] = row
        print(f"{model:22s} " + " | ".join(f"{k}: {v[0][:30]:30s} {v[1] or ''}" for k, v in row.items()), flush=True)
    dest = Path(AUDIT) / "work" / "verify" / "cdotdot" / "fixcheck.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

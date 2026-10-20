"""End-to-end EFT_minkowski_limit verdict (camb.get_background) on A0 vs P_val (A0 + EFTLambdadotdot supplied from EFTLambda
values before compute_positivity_bounds; see minktach_build_patch.py).  EFT_mass_stability = F, EFTCAMB_feedback_level = 2
(the instability messages printed by EFTTestStability are collected).

Usage: python minktach_e2e.py [outjson] [name=tree ...] [--models m1,m2]
"""

import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ROOT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
            "scratchpad/audit")
PAR = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")
OUT = sys.argv[1] if len(sys.argv) > 1 else str(ROOT / "work/verify/minktach/e2e.json")
TREES = {"A0": ROOT / "trees/A0", "P_val": ROOT / "work/verify/minktach/P_val"}
_extra = [x for x in sys.argv[2:] if "=" in x and not x.startswith("--")]
if _extra:
    TREES = {x.split("=", 1)[0]: Path(x.split("=", 1)[1]) for x in _extra}
MODELS = ["5_quint_exponential", "5_quint_inverse_power", "5_quint_quadratic", "5_quint_axion", "5_quint_cubic",
          "5_Kmouflage_1", "5_ADE_1", "5_hdsk_TG_1", "5_Kmimic_1", "5_Horava_eta_1", "5_scg", "2_PEFT_wDE_wCDM_1"]
for _x in sys.argv[2:]:
    if _x.startswith("--models="):
        MODELS = _x.split("=", 1)[1].split(",")

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
    v = "unstable"
sys.stdout.flush()
print("\nRESULT " + json.dumps({"verdict": v, "camb": camb.__file__}), flush=True)
"""


def main() -> None:
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
                                               preexec_fn=lambda: os.nice(10), timeout=900)
                        res = [json.loads(x[7:]) for x in p.stdout.splitlines() if x.startswith("RESULT ")]
                        msgs = sorted({x.strip() for x in p.stdout.splitlines()
                                       if re.search(r"instability|Minkowski|Positivity|Nan|criteria", x)})
                        r = {"model": m, "stab_time": stime, "minkowski": mink, "tree": tname,
                             "verdict": res[0]["verdict"] if res else "error rc=%d" % p.returncode,
                             "camb": res[0]["camb"] if res else None, "messages": msgs}
                        assert r["camb"] is None or str(tree) in r["camb"], r["camb"]
                        rows.append(r)
                        print(f"{m:24s} stab_time={stime:8s} minkowski={mink} {tname:6s}: {r['verdict']:9s} {msgs}",
                              flush=True)
    Path(OUT).write_text(json.dumps(rows, indent=1))
    print("->", OUT)


if __name__ == "__main__":
    main()

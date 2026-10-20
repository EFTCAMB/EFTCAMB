"""Check the proposed fix of eft-pickle-repr-loses-model with EFT-sensitive observables.

For each model, in a subprocess on the chosen tree: reference get_results vs pickle / deepcopy / repr-eval round trip
(sigma8, full lensed TT/EE/TE Cl array up to lmax, H(z)). Also prints the EFT effect size (ref vs EFTflag=0) so a
"restored but reverted to GR" outcome would show up as a large difference.

Usage: python verify_pickle_fix_obs.py [--tree PATH] [--models a,b,c]
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
REPO = Path(__file__).resolve().parents[5]
PARAMS = REPO / "fortran/eftcamb_test/parameters"
WD = Path(f"{ROOT}/work/verify/pickle_fix_obs")

CHILD = r"""
import sys, json, pickle, copy
import numpy as np
tree, ini = sys.argv[1:3]
sys.path.insert(0, tree)
import camb

def obs(p):
    r = camb.get_results(p)
    cl = r.get_cmb_power_spectra(p, CMB_unit="muK", raw_cl=True)["total"][2:2001]
    return {"sigma8": float(r.get_sigma8_0()), "H": r.hubble_parameter([0.5, 1.0]).tolist(), "cl": cl.tolist()}

p = camb.read_ini(ini, no_validate=True)
p.WantTransfer = True
p.set_matter_power(redshifts=[0.0], kmax=2.0)
out = {"ref": obs(p)}
for name, f in (("pickle", lambda q: pickle.loads(pickle.dumps(q))), ("deepcopy", copy.deepcopy),
                ("repr_eval", lambda q: eval(repr(q), {"camb": camb}))):
    out[name] = obs(f(p))
g = p.copy()
g.EFTCAMB.EFTflag = 0
out["gr"] = obs(g)
print("JSON " + json.dumps(out), flush=True)
"""


def cmp(a, b):
    import numpy as np

    ca, cb = np.array(a["cl"]), np.array(b["cl"])
    tt = np.max(np.abs(ca[:, 0] - cb[:, 0]) / np.abs(ca[:, 0]))
    ee = np.max(np.abs(ca[:, 1] - cb[:, 1]) / np.abs(ca[:, 1]))
    te = np.max(np.abs(ca[:, 3] - cb[:, 3]) / np.sqrt(ca[:, 0] * ca[:, 1]))
    s8 = abs(a["sigma8"] - b["sigma8"]) / a["sigma8"]
    return {"sigma8": s8, "TT": tt, "EE": ee, "TE": te}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=f"{ROOT}/work/t0-layout/fix")
    ap.add_argument("--models", default="2_PEFT_wDE_wCDM_1,4_Des_mc5e_CPL_1,5_hdsk_TG_1,4_Des_F_R_wCDM_1")
    a = ap.parse_args()
    WD.mkdir(parents=True, exist_ok=True)
    child = WD / "child.py"
    child.write_text(CHILD)
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    summary = {}
    for m in a.models.split(","):
        src = PARAMS / f"{m}.ini"
        if not src.exists():
            print(f"{m}: no ini")
            continue
        ini = WD / f"{m}.ini"
        ini.write_text(f"DEFAULT({src})\noutput_root = {WD}/{m}\nEFTCAMB_write_background = F\nfeedback_level = 0\n"
                       "EFTCAMB_feedback_level = 0\n")
        p = subprocess.run([sys.executable, str(child), a.tree, str(ini)], capture_output=True, text=True, env=env,
                           cwd=WD, timeout=900)
        res = None
        for ln in p.stdout.splitlines():
            if ln.startswith("JSON "):
                res = json.loads(ln[5:])
        if p.returncode or res is None:
            print(f"{m:22s} rc={p.returncode} {p.stderr[-200:]}")
            summary[m] = {"rc": p.returncode}
            continue
        row = {k: cmp(res["ref"], res[k]) for k in ("pickle", "deepcopy", "repr_eval", "gr")}
        summary[m] = row
        for k, v in row.items():
            print(f"{m:22s} ref vs {k:9s} " + " ".join(f"{kk}={vv:.3e}" for kk, vv in v.items()), flush=True)
    (WD / f"summary_{Path(a.tree).name}.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

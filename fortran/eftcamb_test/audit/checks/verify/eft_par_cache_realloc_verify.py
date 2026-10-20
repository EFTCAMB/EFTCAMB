"""Independent verification of finding eft-par-cache-realloc-crash.

Claim: re-using the CAMBparams held by a CAMBdata (results.Params / its copy / results.get_BAO(z) default) aborts the
process at fortran/results.f90:580 ("allocate of already allocated eft_par_cache%grhormass").

This check differs from the reporter's repro:
  * several models (pure EFT, designer, full mapping) and a same-model GR control (EFTflag=0 override, same cosmology),
  * different API paths: results.calc_background(results.Params.copy()) on the SAME CAMBdata,
    camb.get_results(results.Params) (full pipeline), background-only path (no_thermo), plus controls that re-use the
    USER's CAMBparams (must work: the user's eft_par_cache arrays are never allocated),
  * the -g -fcheck=all AD tree on the background-only path (avoids F-002) to get a backtrace naming the line,
  * B1 (pre-rebase) to show whether the latent bug pre-exists and which entry points reach it there.

Usage: python eft_par_cache_realloc_verify.py --tree A0 [--models ...] [--modes ...] [--no-thermo]
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
            "scratchpad/audit")
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
CHECKS = Path(__file__).resolve().parents[1]

CHILD = r'''
import sys, json
tree, mode, ini, nt, gr = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] == "1", sys.argv[5] == "1"
sys.path.insert(0, tree)
sys.path.insert(1, sys.argv[6])
import camb, harness
z = [0.5, 1.0]
cp = camb.read_ini(ini, no_validate=True)
keys = harness.parse_ini(ini)
if hasattr(cp, "EFTCAMB") and cp.EFTCAMB.EFTflag == 0 and int(keys.get("EFTflag", "0")) != 0:
    cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)   # B1: read_ini ignores EFT keys
print("EFTFLAG", int(cp.EFTCAMB.EFTflag), flush=True)
def bg(p, obj=None):
    if obj is None:
        return camb.get_background(p, no_thermo=nt)
    (obj.calc_background_no_thermo if nt else obj.calc_background)(p)
    return obj
res = bg(cp)
if mode == "fresh":
    out = res.hubble_parameter(z)
elif mode == "user_params_twice":          # control: user's own CAMBparams reused
    out = bg(cp).hubble_parameter(z)
elif mode == "get_BAO_user_params":        # control: get_BAO with explicit user params
    out = res.get_BAO(z, cp)[:, 1]
elif mode == "get_BAO_default":            # CAMB v2 default -> self.Params.copy()
    out = res.get_BAO(z)[:, 1]
elif mode == "same_obj_Params_copy":       # same CAMBdata, its own Params copied
    out = bg(res.Params.copy(), res).hubble_parameter(z)
elif mode == "new_obj_Params":             # new CAMBdata, results.Params passed directly
    out = bg(res.Params).hubble_parameter(z)
elif mode == "get_results_Params":         # full pipeline on results.Params.copy()
    r2 = camb.get_results(res.Params.copy())
    out = r2.hubble_parameter(z)
print("RESULT", json.dumps([float(x) for x in out]), flush=True)
'''

MODES = ["fresh", "user_params_twice", "get_BAO_user_params", "get_BAO_default", "same_obj_Params_copy",
         "new_obj_Params", "get_results_Params"]


def run(tree, mode, ini, nt, gr):
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", CHILD, str(tree), mode, ini, "1" if nt else "0", "1" if gr else "0",
                        str(CHECKS)], capture_output=True, text=True, env=env, cwd=os.path.dirname(ini), timeout=900)
    text = p.stdout + p.stderr
    res = [ln[7:] for ln in p.stdout.splitlines() if ln.startswith("RESULT")]
    flag = [ln.split()[1] for ln in p.stdout.splitlines() if ln.startswith("EFTFLAG")]
    err = [ln.strip() for ln in text.splitlines() if "Fortran runtime error" in ln or ln.startswith("At line")
           or "results.f90" in ln or "Error termination" in ln]
    return {"rc": p.returncode, "H": json.loads(res[0]) if res else None, "EFTflag": flag[0] if flag else None,
            "err": err[:12], "tail": text.splitlines()[-25:] if p.returncode else []}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default="A0")
    ap.add_argument("--models", nargs="+", default=["4_Des_F_R_LCDM_1", "2_PEFT_Omega_lin_1", "5_Kmouflage_1"])
    ap.add_argument("--modes", nargs="+", default=MODES)
    ap.add_argument("--no-thermo", action="store_true")
    ap.add_argument("--gr-control", action="store_true", help="also run each model with EFTflag=0 override")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    tree = Path(a.tree) if os.path.isabs(a.tree) else ROOT / "trees" / a.tree
    tmp = tempfile.mkdtemp(prefix="realloc_verify_", dir=ROOT / "work" / "verify")
    report = {"tree": str(tree), "no_thermo": a.no_thermo, "runs": {}}
    for model in a.models:
        for gr in ([False, True] if a.gr_control else [False]):
            tag = f"{model}{'_GR' if gr else ''}"
            ini = os.path.join(tmp, f"{tag}.ini")
            extra = "EFTflag = 0\n" if gr else ""
            Path(ini).write_text(f"DEFAULT({PARAMS / (model + '.ini')})\noutput_root = {tmp}/{tag}\n"
                                 f"EFTCAMB_write_background = F\nfeedback_level = 0\n{extra}")
            for mode in a.modes:
                r = run(tree, mode, ini, a.no_thermo, gr)
                report["runs"][f"{tag}:{mode}"] = r
                status = "OK " if r["rc"] == 0 else f"rc={r['rc']}"
                msg = r["H"] if r["rc"] == 0 else " | ".join(r["err"][:2])
                print(f"{tree.name:4s} {tag:24s} EFTflag={r['EFTflag']} {mode:22s} {status} {msg}", flush=True)
    if a.out:
        Path(a.out).write_text(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()

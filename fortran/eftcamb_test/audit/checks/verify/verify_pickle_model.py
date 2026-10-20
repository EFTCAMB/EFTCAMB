"""Adversarial verification of finding eft-pickle-repr-loses-model.

Independent of the reporter's repro: builds the EFT model through the python API (camb.set_params with EFT keys)
as well as through read_ini, and checks which serialisation routes keep the EFTCAMB model:
reference, cp.copy(), pickle, copy.deepcopy, repr/eval, multiprocessing (ProcessPoolExecutor, spawn), and a
GR (EFTflag=0) pickle control. Each case runs in its own subprocess with PYTHONPATH=<tree>.

Usage: python verify_pickle_model.py [--tree PATH] [--out JSON]
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
PARAMS = Path(__file__).resolve().parents[3] / "parameters"

PEFT = dict(
    H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965, lmax=1500,
    EFTflag=1, PureEFTmodel=1, PureEFTHorndeski=False, EFTwDE=0, PureEFTmodelOmega=1,
    PureEFTmodelGamma1=0, PureEFTmodelGamma2=0, PureEFTmodelGamma3=0, PureEFTmodelGamma4=0,
    PureEFTmodelGamma5=0, PureEFTmodelGamma6=0, EFTOmega0=0.1,
    EFT_ghost_stability=True, EFT_gradient_stability=True, EFT_mass_stability=False,
    EFT_additional_priors=True,
)

CHILD = r'''
import sys, pickle, copy, json, os
import camb
import numpy as np

def build(src, ini, peft):
    if src == "ini":
        return camb.read_ini(ini, no_validate=True)
    if src == "py":
        return camb.set_params(**peft)
    if src == "gr":
        p = {k: v for k, v in peft.items() if not k.startswith(("EFT", "PureEFT"))}
        p["EFTflag"] = 0
        return camb.set_params(**p)

def observe(p):
    r = camb.get_results(p)
    tt = r.get_cmb_power_spectra(p, CMB_unit="muK")["total"][:, 0]
    return {"H": r.hubble_parameter([0.5, 1.0]).tolist(), "TT10": float(tt[10]), "TT1000": float(tt[1000])}

def worker(p):
    return observe(p)

if __name__ == "__main__":
    mode, src, ini, peft = sys.argv[1], sys.argv[2], sys.argv[3], json.loads(sys.argv[4])
    cp = build(src, ini, peft)
    if mode == "reference":
        p2 = cp
    elif mode == "copy":
        p2 = cp.copy()
    elif mode == "pickle":
        p2 = pickle.loads(pickle.dumps(cp))
    elif mode == "deepcopy":
        p2 = copy.deepcopy(cp)
    elif mode == "repr_eval":
        p2 = eval(repr(cp), {"camb": camb})
    elif mode == "state":
        st = cp.EFTCAMB.__getstate__()
        print("STATE_KEYS", len(st), sorted(st)[:3], "has_model_info", any("model" in k.lower() and k != "EFTCAMB_model_is_designer" for k in st), flush=True)
        print("REPR_HAS_MODEL_NAME", "Pure EFT" in repr(cp) or "f(R)" in repr(cp), flush=True)
        sys.exit(0)
    elif mode == "pool":
        from concurrent.futures import ProcessPoolExecutor
        import multiprocessing as mp
        try:
            with ProcessPoolExecutor(1, mp_context=mp.get_context("spawn")) as ex:
                res = ex.submit(worker, cp).result(timeout=300)
            print("RESULT", json.dumps(res), flush=True)
        except Exception as e:
            print("RESULT exception", type(e).__name__, str(e)[:160], flush=True)
        sys.exit(0)
    print("STATE EFTflag=%d" % p2.EFTCAMB.EFTflag, flush=True)
    try:
        print("RESULT", json.dumps(observe(p2)), flush=True)
    except Exception as e:
        print("RESULT exception", type(e).__name__, str(e)[:160], flush=True)
'''


def run(tree, mode, src, ini, childfile):
    env = dict(os.environ, OMP_NUM_THREADS="2", PYTHONPATH=tree)
    p = subprocess.run([sys.executable, childfile, mode, src, ini, json.dumps(PEFT)], capture_output=True, text=True,
                       env=env, cwd=os.path.dirname(ini), timeout=900)
    lines = [ln for ln in p.stdout.splitlines() if ln.startswith(("STATE", "RESULT", "REPR"))]
    return p.returncode, lines, p.stderr[-600:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=f"{ROOT}/trees/A0")
    ap.add_argument("--out", default=f"{ROOT}/work/verify/pickle_model/verify_A0.json")
    ap.add_argument("--modes", default="reference,copy,pickle,deepcopy,repr_eval,pool,state")
    ap.add_argument("--srcs", default="py,ini,gr")
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix="vpickle_", dir=f"{ROOT}/work/verify/pickle_model")
    childfile = os.path.join(tmp, "child_pickle.py")
    Path(childfile).write_text(CHILD)
    ini = os.path.join(tmp, "fR.ini")
    Path(ini).write_text(f"DEFAULT({PARAMS / '4_Des_F_R_LCDM_1.ini'})\noutput_root = {tmp}/fR\n"
                         "EFTCAMB_write_background = F\nfeedback_level = 0\n")
    out = {}
    print(f"tree: {a.tree}")
    for src in a.srcs.split(","):
        for mode in a.modes.split(","):
            rc, lines, err = run(a.tree, mode, src, ini, childfile)
            status = "OK" if rc == 0 else (f"CRASH signal {-rc}" if rc < 0 else f"CRASH rc={rc}")
            out[f"{src}:{mode}"] = {"rc": rc, "lines": lines, "stderr_tail": err if rc != 0 else ""}
            print(f"{src:3s} {mode:10s} {status:16s} {' | '.join(lines)}", flush=True)
            if rc != 0 and err.strip():
                print("    stderr:", err.strip().splitlines()[-3:])
    Path(a.out).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

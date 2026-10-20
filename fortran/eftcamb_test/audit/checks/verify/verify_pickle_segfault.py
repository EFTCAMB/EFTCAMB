"""Independent adversarial check of finding pickle-deepcopy-eftcamb-segfault.

Differs from the reporter's repro: the EFT models come from read_ini (designer f(R), full-mapping quintessence)
and from EFTCAMB.initialize_parameters with a dict, the crash is probed through several entry points
(EFTCAMB.num_params -> 11_EFTCAMB_python.f90 getter; camb.get_results -> CAMBdata_SetParams; a spawn
multiprocessing worker receiving the params by pickle; repr/eval), the Fortran allocation status of `model` is
read through the guarded getter (model_name returns '' only when .not. allocated(model)), and controls are run:
EFTflag=0 pickle, and a plain upstream (B0) DarkEnergyPPF pickle.

Every case runs in its own subprocess holding one machine-wide CPU slot (harness.slot), OMP_NUM_THREADS=1.

Usage: python verify_pickle_segfault.py --tree A0|B1|B0|<abs path> [--out JSON]
"""

from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
WORK = ROOT / "work" / "verify" / "pickle_segfault"
PARAMS = HERE.parents[2] / "parameters"

FR_DICT = {"EFTflag": 3, "DesignerEFTmodel": 1, "EFTwDE": 0, "EFTB0": 1.0, "EFT_ghost_stability": True,
           "EFT_gradient_stability": True, "EFTCAMB_feedback_level": 0}

CHILD = r'''
import sys, os, json, copy, pickle
tree, source, route, action, payload = sys.argv[1:6]
sys.path.insert(0, tree)
import camb
import numpy as np
payload = json.loads(payload)

XK = {"YHe": payload["yhe"]} if "yhe" in payload else {}

def build():
    if source == "ini":
        cp = camb.read_ini(payload["ini"], no_validate=True)
        return cp
    if source == "dict":
        cp = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965, lmax=1200, **XK)
        cp.EFTCAMB.initialize_parameters(cp, payload["eft"], print_header=False)
        return cp
    if source == "gr":
        return camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965, lmax=1200, **XK)
    if source == "ppf":
        cp = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965, lmax=1200,
                             dark_energy_model="ppf", w=-0.9, wa=0.1, **XK)
        return cp
    raise SystemExit("bad source")

def eft_state(c):
    if not hasattr(c, "EFTCAMB"):
        return {}
    e = c.EFTCAMB
    try:
        name = e.model_name()
    except Exception as ex:  # noqa
        name = "<err %s>" % ex
    return {"EFTflag": int(e.EFTflag), "model_name": name}

cp = build()
if route == "orig":
    c2 = cp
elif route == "copy":
    c2 = cp.copy()
elif route == "deepcopy":
    c2 = copy.deepcopy(cp)
elif route == "pickle":
    c2 = pickle.loads(pickle.dumps(cp))
elif route == "repr":
    import camb.model, camb.dark_energy, camb.initialpower, camb.nonlinear, camb.reionization, camb.recombination
    import camb.sources
    try:
        import camb.eftcamb
    except ImportError:
        pass
    c2 = eval(repr(cp), {"camb": camb})
elif route == "eftonly_deepcopy":
    # deepcopy only the EFTCAMB sub-object and assign it back into a Fortran-level copy
    c2 = cp.copy()
    c2.EFTCAMB = copy.deepcopy(cp.EFTCAMB)
elif route == "mp":
    c2 = cp
else:
    raise SystemExit("bad route")
print("STATE " + json.dumps(eft_state(c2)), flush=True)

out = {}
if action == "num_params":
    out["num_params"] = c2.EFTCAMB.num_params()
elif action == "background":
    r = camb.get_background(c2)
    out["H1"] = float(r.hubble_parameter(1.0))
elif action == "results":
    r = camb.get_results(c2)
    cl = r.get_total_cls(1000)
    out["TT1000"] = float(cl[1000, 0])
    out["sigma8"] = float(r.get_sigma8_0())
elif action == "mp_background":
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    sys.path.insert(0, payload["workerdir"])
    import _vps_worker
    import __main__
    del __main__.__file__  # keep spawn from re-running this script in the worker
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(1, mp_context=ctx, initializer=_vps_worker.init, initargs=(tree,)) as ex:
        try:
            out["H1"] = ex.submit(_vps_worker.background, c2).result(timeout=300)
        except Exception as e:  # BrokenProcessPool when the worker dies
            out["worker_error"] = type(e).__name__ + ": " + str(e)[:200]
print("RESULT " + json.dumps(out), flush=True)
'''

WORKER = '''
import sys

def init(tree):
    sys.path.insert(0, tree)

def background(cp):
    import camb
    return float(camb.get_background(cp).hubble_parameter(1.0))
'''


def make_ini(model: str) -> str:
    ini = WORK / f"{model}.ini"
    ini.write_text(
        f"DEFAULT({PARAMS / (model + '.ini')})\nEFTCAMB_write_background = F\nfeedback_level = 0\n"
        f"EFTCAMB_feedback_level = 0\noutput_root = {WORK / model}\n"
    )
    return str(ini)


def run_case(tree: str, source: str, route: str, action: str, payload: dict) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run(
            [sys.executable, str(WORK / "_vps_child.py"), tree, source, route, action, json.dumps(payload)],
            env=env, capture_output=True, text=True, cwd=str(WORK), timeout=600,
            preexec_fn=lambda: os.nice(10),
        )
    rec = {"source": source, "route": route, "action": action, "label": payload.get("label", "")}
    for line in p.stdout.splitlines():
        if line.startswith("STATE "):
            rec["state"] = json.loads(line[6:])
        elif line.startswith("RESULT "):
            rec["result"] = json.loads(line[7:])
    if p.returncode < 0:
        rec["status"] = "killed:" + signal.Signals(-p.returncode).name
    elif p.returncode > 0:
        rec["status"] = f"rc={p.returncode}"
        rec["stderr_tail"] = p.stderr.strip()[-300:]
    else:
        rec["status"] = "ok"
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default="A0")
    ap.add_argument("--out", default=None)
    ap.add_argument("--quick", action="store_true", help="skip get_results cases")
    args = ap.parse_args()
    tree = args.tree if os.path.isabs(args.tree) else str(ROOT / "trees" / args.tree)
    tag = Path(tree).name
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "_vps_child.py").write_text(CHILD)
    (WORK / "_vps_worker.py").write_text(WORKER)

    sources = []
    if tag != "B0":
        sources.append(("dict", {"eft": FR_DICT, "label": "fR dict"}))
        if tag != "B1":  # B1's read_ini ignores EFTCAMB keys (known fact)
            sources.append(("ini", {"ini": make_ini("4_Des_F_R_LCDM_1"), "label": "4_Des_F_R_LCDM_1"}))
            sources.append(("ini", {"ini": make_ini("5_quint_exponential"), "label": "5_quint_exponential"}))
        sources.append(("gr", {"label": "EFTflag=0 control"}))
    else:
        sources.append(("gr", {"label": "B0 GR control"}))
        sources.append(("ppf", {"label": "B0 PPF control"}))

    recs = []
    for source, payload in sources:
        payload["workerdir"] = str(WORK)
        if tag == "B1":
            payload["yhe"] = 0.2454  # B1 (CAMB 1.6.5) BBN YHe fails with this numpy
        cases = [("orig", "background"), ("copy", "background"), ("deepcopy", "background"),
                 ("pickle", "background"), ("repr", "background"), ("mp", "mp_background")]
        if source in ("dict", "ini"):
            cases += [("pickle", "num_params"), ("orig", "num_params"), ("eftonly_deepcopy", "background")]
            if not args.quick and payload.get("label") == "4_Des_F_R_LCDM_1":
                cases += [("orig", "results"), ("pickle", "results")]
        for route, action in cases:
            rec = run_case(tree, source, route, action, payload)
            recs.append(rec)
            print(json.dumps(rec), flush=True)
    out = args.out or str(WORK / f"verify_{tag}.json")
    Path(out).write_text(json.dumps(recs, indent=1))
    print("wrote", out)


if __name__ == "__main__":
    main()

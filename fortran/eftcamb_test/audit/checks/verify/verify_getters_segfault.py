"""Independent adversarial check of finding eft-getters-gr-segfault.

Differs from the reporter's repro:
  * GR objects are produced by several routes (bare CAMBparams, camb.set_params with no EFT keys, which calls
    initialize_parameters with EFTflag=0, and an in-place re-initialisation from a working f(R) model to GR);
  * positive control: the same getters on an allocated designer f(R) model must return values (rules out a
    ctypes pointer/layout problem as the cause);
  * every failure branch of allocate_EFTCAMB_model (bad EFTflag / Pure / AltPar / Designer / FullMapping model
    index) plus a failure *after* model allocation (bad PureEFTmodelOmega) is probed, recording EFTflag, the
    allocation status of `model` (model_name() returns '' only when .not. allocated(model)), num_params and
    get_results in separate processes;
  * contrast route: camb.read_ini with the same invalid model (camb.f90 path) must raise a clean python error;
  * attribution: lldb backtraces on the -g -O0 -fcheck=all tree (AD) give the faulting Fortran file:line.

Every case runs in its own subprocess holding one machine-wide CPU slot (harness.slot), OMP_NUM_THREADS=1.

Usage: python verify_getters_segfault.py [--tree A0] [--no-lldb] [--out JSON]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import signal
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
WORK = ROOT / "work" / "verify" / "getters_segv"
PARAMS = HERE.parents[2] / "parameters"
PY = sys.executable

FR = {"EFTflag": 3, "DesignerEFTmodel": 1, "EFTwDE": 0, "EFTB0": 1.0, "EFTCAMB_feedback_level": 0}
BAD = {
    "bad_EFTflag": {"EFTflag": 7},
    "bad_Pure": {"EFTflag": 1, "PureEFTmodel": 9},
    "bad_AltPar": {"EFTflag": 2, "AltParEFTmodel": 9},
    "bad_Designer": {"EFTflag": 3, "DesignerEFTmodel": 9},
    "bad_FullMap": {"EFTflag": 4, "FullMappingEFTmodel": 99},
    "bad_Omega_after_alloc": {"EFTflag": 1, "PureEFTmodel": 1, "PureEFTmodelOmega": 99},
}

CHILD = r'''
import sys, json
tree, source, action, payload, debug = sys.argv[1:6]
sys.path.insert(0, tree)
import camb
from ctypes import c_int, byref
payload = json.loads(payload)
debug = debug == "1"

def emit(**kw):
    print("RESULT " + json.dumps(kw), flush=True)

def state(p):
    d = {"EFTflag": int(p.EFTCAMB.EFTflag)}
    if not debug:  # F-001: model_name is rejected by -fcheck on AD
        d["model_name"] = p.EFTCAMB.model_name()
    return d

if source == "bare":
    p = camb.CAMBparams()
elif source == "set_params_gr":
    p = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.12, As=2e-9, ns=0.965)
elif source == "fr":
    p = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.12, As=2e-9, ns=0.965)
    p.EFTCAMB.initialize_parameters(p, payload["fr"], print_header=False)
elif source == "fr_then_gr":
    p = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.12, As=2e-9, ns=0.965)
    p.EFTCAMB.initialize_parameters(p, payload["fr"], print_header=False)
    emit(stage="after_fr", **state(p))
    p.EFTCAMB.initialize_parameters(p, {"EFTflag": 0}, print_header=False)
elif source == "failed":
    if debug:  # F-001: camb.set_params calls read_parameters (string getter rejected by -fcheck on AD)
        p = camb.CAMBparams()
        p.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12)
        p.InitPower.set_params(As=2e-9, ns=0.965)
    else:
        p = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.12, As=2e-9, ns=0.965)
    try:
        p.EFTCAMB.initialize_parameters(p, payload["bad"], print_header=False)
        emit(stage="init", raised=None)
    except Exception as e:
        emit(stage="init", raised=type(e).__name__ + ": " + str(e))
elif source == "read_ini_bad":
    try:
        p = camb.read_ini(payload["ini"], no_validate=True)
        emit(stage="read_ini", raised=None, **state(p))
    except Exception as e:
        emit(stage="read_ini", raised=type(e).__name__ + ": " + str(e)[:120])
        sys.exit(0)
emit(stage="state", **state(p))
p.set_for_lmax(300, lens_potential_accuracy=0)
if action == "none":
    pass
elif action == "num_params_raw":
    # direct ctypes call into 11_EFTCAMB_python.f90, bypassing camb.eftcamb.EFTCAMB.num_params
    n = c_int(-1)
    camb.eftcamb.EFTCAMB_get_num_params(byref(p.EFTCAMB), byref(n))
    emit(stage="action", value=n.value)
elif action in ("num_params", "param_names", "param_values", "param_labels"):
    emit(stage="action", value=getattr(p.EFTCAMB, action)())
elif action == "get_results":
    try:
        r = camb.get_results(p)
        emit(stage="action", value="results", tau0=float(r.tau0),
             TT1000=float(r.get_cmb_power_spectra(CMB_unit="muK")["total"][200, 0]))
    except Exception as e:
        emit(stage="action", value="exception", exc=type(e).__name__ + ": " + str(e)[:150])
'''


def run_child(tree: Path, source: str, action: str, payload: dict, debug: bool = False, timeout: int = 600) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    cmd = [PY, "-c", CHILD, str(tree), source, action, json.dumps(payload), "1" if debug else "0"]
    with harness.slot(1):
        r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=timeout, cwd=str(WORK))
    rc = r.returncode
    res = [json.loads(ln[7:]) for ln in r.stdout.splitlines() if ln.startswith("RESULT ")]
    return {
        "source": source,
        "action": action,
        "rc": rc,
        "signal": signal.Signals(-rc).name if rc < 0 else None,
        "results": res,
        "stderr_tail": r.stderr.strip().splitlines()[-4:],
    }


def lldb_bt(tree: Path, source: str, action: str, payload: dict) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    script = WORK / "lldb_child.py"
    script.write_text(CHILD)
    cmd = ["lldb", "--batch", "-o", "run", "-k", "bt 12", "-k", "quit 1", "--", PY, str(script), str(tree),
           source, action, json.dumps(payload), "1"]
    with harness.slot(1):
        r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=600, cwd=str(WORK))
    out = r.stdout
    stop = [ln.strip() for ln in out.splitlines() if "stop reason" in ln][:1]
    frames = [ln.strip() for ln in out.splitlines() if re.match(r"\s*\*?\s*frame #\d+", ln)]
    return {"source": source, "action": action, "stop": stop, "frames": frames[:12]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default="A0")
    ap.add_argument("--no-lldb", action="store_true")
    ap.add_argument("--out", default=str(WORK / "verify_getters_segfault.json"))
    a = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    tree = harness.tree_path(a.tree)
    out: dict = {"tree": str(tree), "cases": [], "lldb": []}

    def rec(source, action, payload=None, debug=False, t=tree):
        c = run_child(t, source, action, payload or {}, debug=debug)
        out["cases"].append(c)
        last = c["results"][-1] if c["results"] else {}
        how = f"SIG {c['signal']}" if c["signal"] else f"rc={c['rc']}"
        print(f"{Path(t).name:4s} {source:22s} {action:15s} {how:14s} " + " | ".join(json.dumps(x) for x in c["results"]),
              flush=True)
        return c, last

    # 1) GR objects, three routes
    for src in ("bare", "set_params_gr", "fr_then_gr"):
        for act in ("none", "num_params_raw", "num_params", "param_names", "param_values", "get_results"):
            rec(src, act, {"fr": FR})
    # 2) positive control: allocated model
    for act in ("num_params_raw", "num_params", "param_names", "param_labels", "param_values"):
        rec("fr", act, {"fr": FR})
    # 3) every failure branch, reuse of the object
    for name, bad in BAD.items():
        for act in ("none", "num_params", "get_results"):
            rec("failed", act, {"bad": bad})
    # 4) contrast: read_ini with an invalid model goes through camb.f90 (eft_error handled)
    bad_ini = WORK / "bad_pure_model.ini"
    bad_ini.write_text(f"DEFAULT({PARAMS / 'base_params.ini'})\nEFTflag = 1\nPureEFTmodel = 9\n"
                       "EFTCAMB_write_background = F\noutput_root = \n")
    rec("read_ini_bad", "num_params", {"ini": str(bad_ini)})
    rec("read_ini_bad", "get_results", {"ini": str(bad_ini)})
    # 5) attribution with lldb on the debug tree (AD, -g -O0 -fcheck=all)
    if not a.no_lldb:
        ad = harness.tree_path("AD")
        for src, act, pl in (("bare", "num_params_raw", {}), ("bare", "param_values", {}),
                             ("failed", "get_results", {"bad": BAD["bad_Pure"]}),
                             ("failed", "get_results", {"bad": BAD["bad_Omega_after_alloc"]})):
            # first without lldb: what does -fcheck=all report?
            rec(src, act, pl, debug=True, t=ad)
            b = lldb_bt(ad, src, act, pl)
            out["lldb"].append(b)
            print("  lldb", src, act, b["stop"], flush=True)
            for f in b["frames"][:8]:
                print("     ", f[:230], flush=True)
    Path(a.out).write_text(json.dumps(out, indent=1))
    print("wrote", a.out)


if __name__ == "__main__":
    main()

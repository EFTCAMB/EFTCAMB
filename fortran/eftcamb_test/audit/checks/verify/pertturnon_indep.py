"""Independent check of finding pert-turn-on-nonpositive-silently-gr (verify track).

Two methods different from the t9 repro (which used a fresh CAMBparams + initialize_parameters dict, Omega0=0.01):

Part A ("readback"): Python dict path, background only (camb.get_background -> CAMBdata_SetParams, where the RGR
replacement happens).  Reads back results.Params.EFTCAMB.EFTCAMB_pert_turn_on after the call, i.e. the value that
equations.f90 will actually use as the turn-on.  Run for trees A0 and B1 (pre-rebase), in a subprocess per tree.

Part B ("ini"): Fortran ini path through harness.run_model (camb.read_ini of a test ini), model
2_PEFT_Omega_const_1 (Omega0 = 1.0, a large EFT effect), full spectra at the test-ini settings plus the EFT
perturbation evolution at k = 0.01, 0.1 /Mpc.  Compares a_on = 0, -0.1, nan, 2 with: default (0.01), GR (EFTflag=0),
explicit "after today" a_on = 1.1, and the skip_RGR variants (0, nan, 2).

Usage: python pertturnon_indep.py readback|ini|report
Outputs under $AUDIT_ROOT/work/verify/pertturnon/.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "pertturnon"
OUT.mkdir(parents=True, exist_ok=True)
INI = harness.PARAMS_DIR / "2_PEFT_Omega_const_1.ini"

READBACK_CHILD = r"""
import json, math, os, sys
sys.path.insert(0, sys.argv[1])
import camb
base = {"EFTflag": 1, "PureEFTmodel": 1, "PureEFTmodelOmega": 1, "EFTOmega0": 0.01, "EFTwDE": 0}
rows = []
for v, skip in ((0.01, False), (1e-3, False), (1e-30, False), (0.0, False), (-0.1, False), (float("nan"), False),
                (2.0, False), (1.0, False), (0.0, True), (float("nan"), True), (2.0, True)):
    p = camb.CAMBparams()
    p.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12)
    d = dict(base, EFTCAMB_turn_on_time=v, EFTCAMB_skip_RGR=skip)
    p.EFTCAMB.initialize_parameters(p, d, print_header=False)
    read_in = p.EFTCAMB.EFTCAMB_pert_turn_on
    try:
        r = camb.get_background(p)
        used = r.Params.EFTCAMB.EFTCAMB_pert_turn_on
        err = ""
    except Exception as e:
        used, err = float("nan"), repr(e)[:80]
    rows.append({"input": repr(v), "skip_RGR": skip, "after_read": repr(read_in), "after_SetParams": repr(used),
                 "error": err})
print("JSON" + json.dumps(rows))
"""


def readback():
    res = {}
    for tree in ("A0", "B1"):
        tp = str(harness.tree_path(tree))
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        with harness.slot(1):
            p = subprocess.run(["nice", "-n", "10", sys.executable, "-c", READBACK_CHILD, tp], env=env,
                               capture_output=True, text=True, timeout=900, cwd=str(OUT))
        line = [s for s in p.stdout.splitlines() if s.startswith("JSON")]
        res[tree] = json.loads(line[-1][4:]) if line else {"rc": p.returncode, "tail": p.stderr[-2000:]}
    (OUT / "readback.json").write_text(json.dumps(res, indent=1))
    for tree, rows in res.items():
        print(f"== {tree}")
        if isinstance(rows, dict):
            print(rows)
            continue
        print(f"{'input':>8s} {'skipRGR':>7s} {'after read':>12s} {'used (after SetParams)':>24s}  error")
        for r in rows:
            print(f"{r['input']:>8s} {str(r['skip_RGR']):>7s} {r['after_read']:>12s} {r['after_SetParams']:>24s}  "
                  f"{r['error']}")


CASES = {
    "default": {},
    "gr": {"EFTflag": "0"},
    "aon_0": {"EFTCAMB_turn_on_time": "0"},
    "aon_m0p1": {"EFTCAMB_turn_on_time": "-0.1"},
    "aon_nan": {"EFTCAMB_turn_on_time": "nan"},
    "aon_2": {"EFTCAMB_turn_on_time": "2.0"},
    "aon_1p1": {"EFTCAMB_turn_on_time": "1.1"},
    "aon_1em3": {"EFTCAMB_turn_on_time": "1.e-3"},
    "skip_aon_0": {"EFTCAMB_turn_on_time": "0", "EFTCAMB_skip_RGR": "T"},
    "skip_aon_nan": {"EFTCAMB_turn_on_time": "nan", "EFTCAMB_skip_RGR": "T"},
    "skip_aon_2": {"EFTCAMB_turn_on_time": "2.0", "EFTCAMB_skip_RGR": "T"},
    "skip_aon_1p1": {"EFTCAMB_turn_on_time": "1.1", "EFTCAMB_skip_RGR": "T"},
}


def ini_runs(names):
    def one(name):
        st = harness.run_model("A0", INI, OUT / f"ini_{name}.npz", overrides=CASES[name], evol_k=[0.01, 0.1],
                               threads=1, timeout=900)
        return name, st["status"], st.get("log_tail", "")[-300:]

    with ThreadPoolExecutor(4) as ex:
        for name, status, tail in ex.map(one, names):
            print(name, status, "" if status == "ok" else tail)


def report():
    runs = {}
    for name in CASES:
        f = OUT / f"ini_{name}.npz"
        js = OUT / f"ini_{name}.json"
        status = json.loads(js.read_text())["status"] if js.exists() else "missing"
        runs[name] = (status, harness.load(f) if f.exists() and status == "ok" else None)
        if runs[name][1] is not None:
            meta = json.loads(str(runs[name][1]["meta"]))
            if meta.get("errors"):
                status += f" errors={meta['errors']}"
                runs[name] = (status, runs[name][1])

    def tt(d):
        return d["cl_total"][2:2500, 0]

    def bb(d):
        return d["cl_tensor"][2:1500, 2] if d["cl_tensor"].shape[0] > 1500 else None

    ref, gr = runs["default"][1], runs["gr"][1]
    print(f"model {INI.name}; GR vs default max|dTT/TT| = {np.max(np.abs(tt(gr) / tt(ref) - 1)):.4e}")
    print(f"{'case':14s} {'status':28s} {'|TT/TTdef-1|':>13s} {'|TT/TTgr-1|':>12s} {'|TT/TT1p1-1|':>13s} "
          f"{'bit=1p1':>8s} {'max|pi| k=.1':>13s} {'tensBB/1p1-1':>13s}")
    one = runs["aon_1p1"][1]
    for name, (status, d) in runs.items():
        if d is None:
            print(f"{name:14s} {status[:28]:28s}")
            continue
        c = tt(d)
        pi = np.nan
        for key in ("evol_pi", "evol_EFTpi", "evol_pi_field"):
            if key in d:
                pi = np.nanmax(np.abs(d[key][1]))
                break
        tb = np.nan
        if bb(d) is not None and bb(one) is not None:
            tb = np.max(np.abs(bb(d) / bb(one) - 1))
        bit = all(np.array_equal(d[k], one[k], equal_nan=True) for k in ("cl_total", "cl_tensor", "pk"))
        print(f"{name:14s} {status[:28]:28s} {np.max(np.abs(c / tt(ref) - 1)):13.4e} "
              f"{np.max(np.abs(c / tt(gr) - 1)):12.4e} {np.max(np.abs(c / tt(one) - 1)):13.4e} {str(bit):>8s} "
              f"{pi:13.4e} {tb:13.4e}")
    d = runs["default"][1]
    print("evol fields:", [k for k in d.keys() if k.startswith("evol_")][:40])


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "readback":
        readback()
    elif mode == "ini":
        ini_runs(sys.argv[2:] or list(CASES))
    elif mode == "report":
        report()

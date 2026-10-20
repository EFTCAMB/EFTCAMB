"""Verify track, finding pert-turn-on-nonpositive-silently-gr: separate the two regimes of EFTCAMBReturnToGR for an
invalid turn-on (a_on = 0, <0, NaN):
  - EFT functions constant in a (Omega const): f(a=NaN) is finite > threshold -> RGR_time = a_used = NaN;
  - EFT functions depending on a (Omega linear):  f(NaN) = NaN -> loop never fires -> RGR_time = 1.1.
For each, run the ini path (harness, A0) with tensors on (test-ini default) and off, and compare with GR (EFTflag=0),
the default turn-on 0.01, and the explicit "after today" a_on = 1.1.  Also reads back the turn-on value used after
CAMBdata_SetParams for both models (Python path, background only).

Usage: python pertturnon_regimes.py run|report
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

OUT = harness.AUDIT_ROOT / "work" / "verify" / "pertturnon" / "regimes"
OUT.mkdir(parents=True, exist_ok=True)
MODELS = {"const": "2_PEFT_Omega_const_1.ini", "lin": "2_PEFT_Omega_lin_1.ini"}
AON = {"def": None, "gr": "GR", "0": "0", "m0p1": "-0.1", "nan": "nan", "1p1": "1.1", "2": "2.0"}
TENS = {"tensT": "T", "tensF": "F"}

CHILD = r"""
import json, sys
sys.path.insert(0, sys.argv[1])
import camb
rows = []
for om in (1, 2):
    for v in (0.01, 0.0, -0.1, float("nan"), 2.0):
        p = camb.CAMBparams()
        p.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12)
        d = {"EFTflag": 1, "PureEFTmodel": 1, "PureEFTmodelOmega": om, "EFTOmega0": 1.0, "EFTOmegaExp": 1.0,
             "EFTwDE": 0, "EFTCAMB_turn_on_time": v}
        p.EFTCAMB.initialize_parameters(p, d, print_header=False)
        r = camb.get_background(p)
        rows.append([om, repr(v), repr(r.Params.EFTCAMB.EFTCAMB_pert_turn_on)])
print("JSON" + json.dumps(rows))
"""


def cases():
    for m in MODELS:
        for t in TENS:
            for a in AON:
                yield m, t, a


def run():
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run(["nice", "-n", "10", sys.executable, "-c", CHILD, str(harness.TREES["A0"])], env=env,
                           capture_output=True, text=True, timeout=600, cwd=str(OUT))
    rows = json.loads([s for s in p.stdout.splitlines() if s.startswith("JSON")][-1][4:])
    (OUT / "readback.json").write_text(json.dumps(rows))
    print("readback (PureEFTmodelOmega, input a_on, a_on used after CAMBdata_SetParams):")
    for r in rows:
        print("  ", r)

    def one(c):
        m, t, a = c
        ov = {"get_tensor_cls": TENS[t]}
        if AON[a] == "GR":
            ov["EFTflag"] = "0"
        elif AON[a] is not None:
            ov["EFTCAMB_turn_on_time"] = AON[a]
        st = harness.run_model("A0", harness.PARAMS_DIR / MODELS[m], OUT / f"{m}_{t}_{a}.npz", overrides=ov,
                               threads=1, timeout=900)
        err = ""
        if st["status"] != "ok":
            tail = st.get("log_tail", "")
            tail = tail if isinstance(tail, str) else "\n".join(tail)
            err = [s for s in tail.splitlines() if "CAMBError" in s or "Error in" in s or "EFTCAMB" in s][-1:]
        return c, st["status"], err

    with ThreadPoolExecutor(4) as ex:
        for c, status, err in ex.map(one, list(cases())):
            print(c, status, err)


def report():
    def ld(m, t, a):
        js = OUT / f"{m}_{t}_{a}.json"
        if not js.exists():
            return "missing", None
        st = json.loads(js.read_text())["status"]
        return st, (harness.load(OUT / f"{m}_{t}_{a}.npz") if st == "ok" else None)

    for m in MODELS:
        for t in TENS:
            _, ref = ld(m, t, "def")
            _, gr = ld(m, t, "gr")
            _, one = ld(m, t, "1p1")
            print(f"== {MODELS[m]} tensors={TENS[t]}: GR vs default max|dTT/TT| (l<2500) = "
                  f"{np.max(np.abs(gr['cl_total'][2:2500, 0] / ref['cl_total'][2:2500, 0] - 1)):.4e}")
            print(f"   {'a_on':6s} {'status':14s} {'|TT/TTdef-1|':>13s} {'|TT/TTgr-1|':>12s} {'|TT/TT1.1-1|':>13s} "
                  f"{'|BBt/BBt1.1-1|':>15s} {'bit==1.1':>8s}")
            for a in AON:
                st, d = ld(m, t, a)
                if d is None:
                    print(f"   {a:6s} {st:14s}")
                    continue
                c = d["cl_total"][2:2500, 0]
                if TENS[t] == "T":
                    bt = np.max(np.abs(d["cl_tensor"][2:1500, 2] / one["cl_tensor"][2:1500, 2] - 1))
                else:
                    bt = np.nan
                bit = all(np.array_equal(d[k], one[k], equal_nan=True) for k in ("cl_total", "cl_tensor", "pk"))
                print(f"   {a:6s} {st:14s} {np.max(np.abs(c / ref['cl_total'][2:2500, 0] - 1)):13.4e} "
                      f"{np.max(np.abs(c / gr['cl_total'][2:2500, 0] - 1)):12.4e} "
                      f"{np.max(np.abs(c / one['cl_total'][2:2500, 0] - 1)):13.4e} {bt:15.4e} {str(bit):>8s}")


if __name__ == "__main__":
    {"run": run, "report": report}[sys.argv[1]]()

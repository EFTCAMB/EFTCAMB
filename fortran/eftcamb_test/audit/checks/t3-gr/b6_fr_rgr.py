"""t3-gr (b): designer f(R) approach to GR across the RGR switch.

fR__s1em06 (B0=1e-6) has RGR time 0.366 and P(k) +55% at k=30/Mpc, B0=1e-8 has RGR time 1.1 (pure GR). Scan B0 in
between with the default EFTCAMB_GR_threshold (1e-8) and with the RGR switch disabled (threshold=-1: EFT from the
turn-on time) to measure the error the RGR approximation introduces near GR, and whether the approach is smooth.
Usage: python b6_fr_rgr.py [default|forced] (runs that chunk; then prints the table of whatever exists)
"""

import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from common import PARAMS, RUNS, harness

INI = PARAMS / "4_Des_F_R_LCDM_1.ini"
SUB = RUNS / "b6"
REF = RUNS / "b1" / "REF_base__A0g.npz"
VALS = [1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-9]
MODES = {"default": {}, "forced": {"EFTCAMB_GR_threshold": "-1"}}


def name(mode, v):
    return f"fR_{mode}_{v:g}".replace("-", "m")


def one(job):
    mode, v = job
    out = SUB / f"{name(mode, v)}.npz"
    if not out.with_suffix(".json").exists():
        harness.run_model("A0", INI, out, overrides={"EFTB0": repr(v), **MODES[mode]}, threads=1, timeout=1500)
    return job


def main() -> None:
    chunk = sys.argv[1:] or []
    jobs = [(m, v) for m in chunk for v in VALS]
    with ThreadPoolExecutor(3) as ex:
        list(ex.map(one, jobs))
    ref = harness.load(REF)
    k = ref["pk_k"]
    ik = [int(np.argmin(abs(k - x))) for x in (0.1, 1.0, 3.0, 10.0, 30.0)]
    print(f"{'mode':8s} {'B0':>7s} {'status':10s} {'RGR a':>8s} " + " ".join(f"P(k={k[i]:.2g})" for i in ik)
          + "   sigma8   max|dTT/TT|")
    for m in MODES:
        for v in VALS:
            out = SUB / f"{name(m, v)}.npz"
            js = out.with_suffix(".json")
            if not js.exists():
                continue
            st = json.loads(js.read_text())["status"]
            if st != "ok":
                tail = [ln for ln in out.with_suffix(".log").read_text(errors="replace").splitlines()
                        if "nstab" in ln or "rror" in ln][-1:]
                print(f"{m:8s} {v:7.0e} {st:10s} {tail}")
                continue
            d = harness.load(out)
            rgr = re.findall(r"Return to GR time: *([-+0-9.E]+)", out.with_suffix(".log").read_text())
            r = d["pk"][0] / ref["pk"][0] - 1
            c = harness.compare(out, REF)
            print(f"{m:8s} {v:7.0e} {st:10s} {float(rgr[-1]) if rgr else float('nan'):8.4f} "
                  + " ".join(f"{r[i]:+11.2e}" for i in ik)
                  + f" {d['sigma8'][0] / ref['sigma8'][0] - 1:+9.1e} {c['cl']['total']['TT']:9.1e}")


if __name__ == "__main__":
    main()

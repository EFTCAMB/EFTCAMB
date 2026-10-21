"""Verifier (spline-integral-convention): can the broken cubic integral() run SILENTLY with a smooth node layout?
30 log-spaced nodes 1e-9..1.2 (cubic spline of yint smooth and positive), constant w = -0.9, designer f(R) and
shift-symmetric alpha_B on A0; compare with the EFTwDE = 1 references of verify_spline_integral_models.py.
usage: /Users/marco/Software/miniforge3/bin/python verify_spline_integral_silent.py
"""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
from verify_spline_integral_models import DFR, OUT, SS, spline  # noqa: E402

NODES = list(np.logspace(-9, np.log10(1.2), 30))
CASES = {"DFR_s9log": (DFR, "DFR_ref"), "SS_s9log": (SS, "SS_ref")}


def run(name):
    return name, harness.run_model("A0", CASES[name][0], OUT / f"{name}.npz", overrides=spline(9, NODES), threads=1,
                                   timeout=900)


with ThreadPoolExecutor(2) as ex:
    for name, st in ex.map(run, CASES):
        msg = ""
        if st["status"] != "ok":
            msg = " / ".join(ln.strip()[:120] for ln in st["log_tail"] if ln.strip())[-400:]
            print(name, st["status"], msg)
            continue
        d = harness.load(OUT / f"{name}.npz")
        r = harness.load(OUT / f"{CASES[name][1]}.npz")
        c = harness.compare(d, r)
        a = r["eft_a_grid"]
        i = np.searchsorted(a, 1e-3)
        print(f"{name} ok: grhov_t(1e-3) test {d['eft_grhov_t'][i]:.4e} ref {r['eft_grhov_t'][i]:.4e}; "
              f"grhov_t(1) test {d['eft_grhov_t'][-1]:.4e} ref {r['eft_grhov_t'][-1]:.4e}; worst Cl {c['worst_cl']:.3e}"
              f" sigma8 {c['sigma8']:.3e}")

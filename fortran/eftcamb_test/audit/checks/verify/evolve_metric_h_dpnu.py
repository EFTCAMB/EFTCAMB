"""Verifier: does using the true massive-neutrino pressure perturbation (patched tree emh_dpnu) remove the
evolve_metric_h low-ell TT discrepancy?"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "verify"
OUT = W / "evolve_metric_h"
TREE = W / "emh_dpnu"
INI = harness.PARAMS_DIR / "5_hdsk_TG_1.ini"
jobs = [(f"dpnu_{v}_{h}", ov | {"EFTCAMB_evolve_metric_h": h}) for v, ov in
        (("full", {}), ("nocubic", {"Horndeskic30": "0"})) for h in ("F", "T")]


def run(j):
    return j[0], harness.run_model(TREE, INI, OUT / f"{j[0]}.npz", overrides=j[1], threads=1)


with ThreadPoolExecutor(4) as ex:
    st = dict(ex.map(run, jobs))
print({n: s["status"] for n, s in st.items()})
ells, ks = [2, 5, 10, 30, 100], [0.1, 1.0, 10.0, 50.0, 100.0]


def cmp(a, b):
    ra, rb = harness.load(OUT / f"{a}.npz"), harness.load(OUT / f"{b}.npz")
    tt = rb["cl_unlensed_scalar"][:200, 0] / ra["cl_unlensed_scalar"][:200, 0] - 1
    pa = np.interp(np.log(ks), np.log(ra["pk_k"]), np.log(ra["pk"][0]))
    pb = np.interp(np.log(ks), np.log(rb["pk_k"]), np.log(rb["pk"][0]))
    print(f"{b} vs {a}: dTT " + " ".join(f"l{e}:{tt[e]:+.2e}" for e in ells) + " | dP " +
          " ".join(f"k{k:g}:{x:+.2e}" for k, x in zip(ks, np.exp(pb - pa) - 1)))


cmp("dpnu_full_F", "dpnu_full_T")
cmp("dpnu_nocubic_F", "dpnu_nocubic_T")
cmp("evol_F", "dpnu_full_F")  # effect of patch on default path

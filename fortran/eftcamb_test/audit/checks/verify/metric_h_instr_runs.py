"""Runs on the instrumented tree (metric_h_build_instr.py) for finding evolve-metric-h-inconsistent (F-065).

Cases (5_hdsk_TG_1, default test settings):
  F        : default constraint path, logged (h'_00 used, h'' from the ii equation evaluated on that trajectory)
  T        : evolve_metric_h=T as shipped (IC h' = GR 2kz), logged
  T_icfix  : evolve_metric_h=T with IC h' from the EFT 00 constraint (EFT_MH_ICFIX=1), logged
  F_nolog / T_nolog : no env vars -> must reproduce A0 bit-for-bit (instrumentation inert)
Then prints spectra differences w.r.t. F.

Usage: [MH_TREE=mh_instr_dgp] python metric_h_instr_runs.py [model_ini_stem] [extra ini overrides as key=value ...]
"""

import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

TREE_NAME = os.environ.get("MH_TREE", "mh_instr")  # e.g. mh_instr_dgp
TREE = harness.AUDIT_ROOT / "work" / "verify" / TREE_NAME
OUT = harness.AUDIT_ROOT / "work" / "verify" / "metric_h"
stem = sys.argv[1] if len(sys.argv) > 1 else "5_hdsk_TG_1"
extra = dict(a.split("=", 1) for a in sys.argv[2:])
tag = stem + ("_" + "_".join(f"{k}{v}" for k, v in extra.items()) if extra else "")
if TREE_NAME != "mh_instr":
    tag = f"{TREE_NAME}_{tag}"
INI = harness.PARAMS_DIR / f"{stem}.ini"
if not INI.exists():
    INI = harness.REPO / "fortran" / "eftcamb_test" / "audit" / "inis" / f"{stem}.ini"

CASES = {
    "F": ({"EFTCAMB_evolve_metric_h": "F"}, {"EFT_MH_LOG": str(OUT / f"{tag}_F.bin")}),
    "T": ({"EFTCAMB_evolve_metric_h": "T"}, {"EFT_MH_LOG": str(OUT / f"{tag}_T.bin")}),
    "T_icfix": ({"EFTCAMB_evolve_metric_h": "T"}, {"EFT_MH_LOG": str(OUT / f"{tag}_Ticfix.bin"), "EFT_MH_ICFIX": "1"}),
}
if not extra and stem == "5_hdsk_TG_1" and TREE_NAME == "mh_instr":
    CASES["T_nolog"] = ({"EFTCAMB_evolve_metric_h": "T"}, {})


def run(item):
    name, (ov, env) = item
    npz = OUT / f"instr_{tag}_{name}.npz"
    if npz.exists():
        return name, {"status": "ok", "cached": True}
    return name, harness.run_model(TREE, INI, npz, overrides={**extra, **ov}, extra_env=env, threads=1)


with ThreadPoolExecutor(4) as ex:
    st = dict(ex.map(run, CASES.items()))
for n, s in st.items():
    if s["status"] != "ok":
        raise SystemExit(f"{n}: {s['status']} {s.get('log_tail', [])[-5:]}")

d = {n: harness.load(OUT / f"instr_{tag}_{n}.npz") for n in CASES}
ells = [2, 4, 5, 10, 20, 30, 100, 1000]
ks = [0.01, 0.1, 1.0, 5.0, 10.0, 20.0, 50.0, 100.0]
ref = d["F"]
res = {}
print(f"model {stem} overrides {extra}; reference F (constraint path)")
print(f"{'case':10s} " + " ".join(f"TT l={e:<5d}" for e in ells) + " | " + " ".join(f"P k={k:<6g}" for k in ks))
for n in CASES:
    if n == "F":
        continue
    x = d[n]
    m = min(len(x["cl_unlensed_scalar"]), len(ref["cl_unlensed_scalar"]))
    tt = x["cl_unlensed_scalar"][:m, 0] / ref["cl_unlensed_scalar"][:m, 0] - 1
    kmax = min(x["pk_k"].max(), ref["pk_k"].max())
    pr = np.interp(np.log(ks), np.log(ref["pk_k"]), np.log(ref["pk"][0]))
    px = np.interp(np.log(ks), np.log(x["pk_k"]), np.log(x["pk"][0]))
    pk = np.where(np.array(ks) <= kmax, np.exp(px - pr) - 1, np.nan)
    res[n] = {"TT": {e: float(tt[e]) for e in ells}, "pk": {k: float(p) for k, p in zip(ks, pk)}}
    print(f"{n:10s} " + " ".join(f"{tt[e]:+10.2e}" for e in ells) + " | " + " ".join(f"{p:+9.2e}" for p in pk))
if "T_nolog" in d:
    a0 = harness.load(harness.AUDIT_ROOT / "work" / "t8a-equations" / "repro_metric_h" / "run1.npz")
    same = np.array_equal(a0["cl_unlensed_scalar"], d["T_nolog"]["cl_unlensed_scalar"]) and np.array_equal(
        a0["pk"], d["T_nolog"]["pk"])
    a0F = harness.load(harness.AUDIT_ROOT / "work" / "t8a-equations" / "repro_metric_h" / "run0.npz")
    sameF = np.array_equal(a0F["cl_unlensed_scalar"], d["F"]["cl_unlensed_scalar"])
    print(f"instrumented T without env == A0 T (repro run1) bitwise: {same};  logged F == A0 F (repro run0): {sameF}")
    res["inert"] = {"T_nolog_eq_A0": bool(same), "F_log_eq_A0": bool(sameF)}
(OUT / f"instr_{tag}_summary.json").write_text(json.dumps(res, indent=1))

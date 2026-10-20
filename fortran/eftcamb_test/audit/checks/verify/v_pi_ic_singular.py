"""Independent check of t3-pi-ic-singular-denominator: does the P(k) spike come from the EFT_IC_type=1 initial condition?
Variants at RPHmassP0=1e-6: IC types 1/2/3, shifted turn-on time, dense k sampling."""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

INI = harness.PARAMS_DIR / "3_RPH_MK_const_1.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "pi_ic"
PK = {"get_scalar_cls": "F", "get_tensor_cls": "F", "do_lensing": "F", "transfer_kmax": "10", "num_redshiftwindows": "0"}
DENSE = {"transfer_k_per_logint": "60"}
JOBS = [
    ("GR", GR, {}, {"Accuracy.IntTolBoost": 2.0}),
    ("GR_dense", GR, DENSE, {"Accuracy.IntTolBoost": 2.0}),
    ("ic1", INI, {"RPHmassP0": "1e-06", "EFT_IC_type": "1"}, {}),
    ("ic2", INI, {"RPHmassP0": "1e-06", "EFT_IC_type": "2"}, {}),
    ("ic3", INI, {"RPHmassP0": "1e-06", "EFT_IC_type": "3"}, {}),
    ("ic1_on005", INI, {"RPHmassP0": "1e-06", "EFTCAMB_turn_on_time": "0.005"}, {}),
    ("ic1_on02", INI, {"RPHmassP0": "1e-06", "EFTCAMB_turn_on_time": "0.02"}, {}),
    ("ic1_dense", INI, {"RPHmassP0": "1e-06", **DENSE}, {}),
    ("ic2_dense", INI, {"RPHmassP0": "1e-06", "EFT_IC_type": "2", **DENSE}, {}),
    ("ic1_m1e-8_dense", INI, {"RPHmassP0": "1e-08", **DENSE}, {}),
]


def job(j):
    name, ini, ov, ps = j
    st = harness.run_model("A0", ini, OUT / f"{name}.npz", overrides={**PK, **ov}, python_set=ps, threads=1, timeout=600)
    return name, st


with ThreadPoolExecutor(4) as ex:
    res = dict(ex.map(job, JOBS))
for n, st in res.items():
    if st["status"] != "ok":
        print(n, st["status"], st.get("log_tail"))
d = {n: harness.load(OUT / f"{n}.npz") for n, st in res.items() if st["status"] == "ok"}
for n in d:
    if n.startswith("GR"):
        continue
    ref = d["GR_dense"] if "dense" in n else d["GR"]
    if len(ref["pk_k"]) != len(d[n]["pk_k"]) or not np.allclose(ref["pk_k"], d[n]["pk_k"]):
        print(n, "k grid mismatch")
        continue
    r = np.abs(d[n]["pk"][0] / ref["pk"][0] - 1)
    i = int(np.argmax(r))
    top = np.argsort(r)[-5:][::-1]
    print(f"{n:16s} max {r[i]:.2e} at k={ref['pk_k'][i]:.4f}  median {np.median(r):.1e}  nk={len(r)}  "
          f"top k: {[f'{ref['pk_k'][t]:.4f}:{r[t]:.1e}' for t in top]}")

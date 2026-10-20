"""Independent check of kmimic-pk-single-mode-glitch.

1. Dense transfer k sampling (transfer_k_per_logint=60): is the feature a narrow window located at k* = 0.001/tau_on?
2. Move EFTCAMB_turn_on_time 1e-8 -> 1e-7: k* should move by ~10x (radiation era tau ~ a); does the feature follow?
3. Compute tau_on directly from the background (camb results.conformal_time) to get k*.
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "kmimic_glitch"
OUT.mkdir(parents=True, exist_ok=True)
BASE = {"num_redshiftwindows": "0", "transfer_k_per_logint": "60"}


def dev(d):
    k, p = d["pk_k"], d["pk"][0]
    m = (k > 1e-3) & (k < 1.5)
    lk, lp = np.log(k[m]), np.log(p[m])
    h1, h2 = lk[1:-1] - lk[:-2], lk[2:] - lk[1:-1]
    d2 = lp[1:-1] - (lp[:-2] * h2 + lp[2:] * h1) / (h1 + h2)
    return k[m][1:-1], d2


def job(a):
    tag, model, over = a
    out = OUT / f"{model}_{tag}.npz"
    st = harness.run_model("A0", harness.PARAMS_DIR / f"{model}.ini", out, overrides={**BASE, **over}, threads=1,
                           timeout=1500)
    if st["status"] != "ok":
        return (model, tag), st["status"] + " " + str(st.get("log_tail", ""))[-600:]
    return (model, tag), harness.load(out)


jobs = [("dense", "1_EFT_GR", {}), ("dense", "5_Kmimic_1", {}),
        ("dense_on1e-7", "5_Kmimic_1", {"EFTCAMB_turn_on_time": "1.e-7"}),
        ("dense_on1e-9", "5_Kmimic_1", {"EFTCAMB_turn_on_time": "1.e-9"})]
with ThreadPoolExecutor(4) as ex:
    res = dict(ex.map(job, jobs))

for key, d in res.items():
    if isinstance(d, str):
        print(key, d)
        continue
    k, d2 = dev(d)
    big = np.where(np.abs(d2) > 0.05)[0]
    print(f"{key}: n_k={len(d['pk_k'])} sigma8={d['sigma8'][0]:.5f} max|dev|={np.abs(d2).max():.3e} "
          f"at k={k[np.argmax(np.abs(d2))]:.4f}; points |dev|>0.05: "
          + ", ".join(f"{k[i]:.4f}:{d2[i]:+.2f}" for i in big[:12]))

# P ratio Kmimic/GR around k* on the dense grid
g, km = res[("1_EFT_GR", "dense")], res[("5_Kmimic_1", "dense")]
if not isinstance(g, str) and not isinstance(km, str):
    kk = km["pk_k"]
    r = km["pk"][0] / np.exp(np.interp(np.log(kk), np.log(g["pk_k"]), np.log(g["pk"][0])))
    sel = (kk > 0.15) & (kk < 0.32)
    print("k, P_Kmimic/P_GR (dense, around k*):")
    print("  " + " ".join(f"{a:.4f}:{b:.3f}" for a, b in zip(kk[sel], r[sel])))
    meta = km["meta"]
    print("meta keys:", list(json.loads(str(meta)).keys()) if not isinstance(meta, dict) else list(meta))

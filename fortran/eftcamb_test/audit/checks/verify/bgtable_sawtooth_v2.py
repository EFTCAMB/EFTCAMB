"""Verifier (t3-bg-table-resolution-gr-floor): independent check via dense sampling of EFT background functions.

For each model, evaluate EFTLambda, adotoa (H) on a dense a-grid around a in [0.45, 0.55] with default and 10x table
points; the default-vs-10x difference must show a sawtooth with period = table spacing in ln a, ~zero at table nodes
and peak ~ h^2/8 |f''|/|f| at midpoints (linear interpolation signature). Also print H deviation from an analytic
LCDM-like reference is NOT used (the 10x run is the reference).
"""
import os
import sys
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(harness.TREES["A0"]))
import camb  # noqa: E402

P = harness.PARAMS_DIR
OUT = harness.AUDIT_ROOT / "work" / "verify" / "bgtable_saw"
MODELS = {
    "Kmouflage": (P / "5_Kmouflage_1.ini", {"eps2_0": "-1e-08"}, 1000, np.log(1e-9), 0.0),
    "ADE": (P / "5_ADE_1.ini", {"f_ac": "1e-08"}, 1500, np.log(1e-10), 0.0),
    "quint": (P / "5_quint_exponential.ini", {"EFT_lambda": "0.003"}, 1100, np.log(1e-8), np.log(1.1)),
}


def run(ini, ov, n, a):
    w = OUT / "tmp.ini"
    lines = [f"DEFAULT({ini})", "EFTCAMB_write_background = F", "get_tensor_cls = F", "do_lensing = F"]
    lines += [f"{k} = {v}" for k, v in ov.items()]
    if n:
        lines.append(f"model_background_num_points = {n}")
    w.write_text("\n".join(lines) + "\n")
    cp = camb.read_ini(str(w), no_validate=True)
    res = camb.get_background(cp)
    fields, arr = cp.EFTCAMB.get_eft_functions(res, a)
    H = res.hubble_parameter(1 / a - 1)
    return {f: np.array(arr[f]) for f in ("EFTLambda", "adotoa", "EFTc")}, H


with harness.slot(1):
    for name, (ini, ov, n0, x0, x1) in MODELS.items():
        h = (x1 - x0) / (n0 - 1)
        a = np.exp(np.linspace(np.log(0.45), np.log(0.55), 2001))
        d0, H0 = run(ini, ov, None, a)
        d1, H1 = run(ini, ov, 10 * n0, a)
        x = np.log(a)
        phase = ((x - x0) / h) % 1.0
        rL = d0["EFTLambda"] / d1["EFTLambda"] - 1
        rH = H0 / H1 - 1
        node = (phase < 0.03) | (phase > 0.97)
        mid = np.abs(phase - 0.5) < 0.03
        # predicted linear-interp error for f ~ a^2: h^2/8 * 4 (sign: chord above convex fn)
        print(f"{name:10s} h={h:.4f} pred_mid(Lambda~a^2)={h*h/2:.2e} | dLambda/Lambda node={np.abs(rL[node]).max():.1e}"
              f" mid={np.abs(rL[mid]).mean():.1e} max={np.abs(rL).max():.1e} | dH/H node={np.abs(rH[node]).max():.1e}"
              f" mid={np.abs(rH[mid]).mean():.1e} max={np.abs(rH).max():.1e}")
        np.savez(OUT / f"{name}.npz", a=a, rL=rL, rH=rH, phase=phase)

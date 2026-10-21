"""Physics-lens verification of finding spline-cubic-endpoints-derivatives (cubic spline parametrization, flag 9).

Questions:
  1. Does the end-node null (value(a=1) = null_value) matter for the spectra, or only the zeroed nodal derivatives?
     -> lin_1 (Omega = a) with nodes 0..0.8 and the last node at 1 + 1e-6 (value(a<=1) never null, derivative zeroing
        still at ~1) compared with the reported N=6 [0,1] case.
  2. Is the effect still large at a realistic amplitude?  lin_2 (Omega = 0.1 a), spline N=6 on the default x nodes.
  3. Curved but exactly representable data: power_2 (Omega = 0.1 a^2), spline N=6 on [0,1] (a cubic spline with
     bctype 2 reproduces any cubic exactly, so value is exact; only the derivative construction is tested).
  4. Pure numpy: for a non-polynomial function, is the "separate spline of nodal derivatives" design itself fine once
     the endpoint zeroing is removed (i.e. is the zeroing the operative defect)?
usage: python verify_spline_cubic_physics.py
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "spline_cubic_phys"
P = harness.PARAMS_DIR


def spline_ov(xs, vs):
    ov = {"PureEFTmodelOmega": "9", "EFTOmega_Spline_Pixels": str(len(xs)), "EFTOmega_null_value": "0.0"}
    for i, (x, v) in enumerate(zip(xs, vs)):
        ov[f"EFTOmegax{i + 1}"] = repr(float(x))
        ov[f"EFTOmegav{i + 1}"] = repr(float(v))
    return ov


x6 = np.linspace(0, 1, 6)
x6e = x6.copy()
x6e[-1] = 1.0 + 1e-6
CASES = {
    "lin1_ref": (P / "2_PEFT_Omega_lin_1.ini", {}),
    "lin1_N6": (P / "2_PEFT_Omega_lin_1.ini", spline_ov(x6, x6)),
    "lin1_N6_xend": (P / "2_PEFT_Omega_lin_1.ini", spline_ov(x6e, x6e)),
    "lin2_ref": (P / "2_PEFT_Omega_lin_2.ini", {}),
    "lin2_N6": (P / "2_PEFT_Omega_lin_2.ini", spline_ov(x6, 0.1 * x6)),
    "pow2_ref": (P / "2_PEFT_Omega_power_2.ini", {}),
    "pow2_N6": (P / "2_PEFT_Omega_power_2.ini", spline_ov(x6, 0.1 * x6**2)),
}


def run(name):
    ini, ov = CASES[name]
    out = OUT / f"{name}.npz"
    if out.exists():
        return name, {"status": "ok(cached)"}
    return name, harness.run_model("A0", ini, out, overrides=ov, threads=1)


def numpy_design_check():
    def spline3(x, y):
        cl = np.polyfit(x[:4] - x[0], y[:4], 3)
        cr = np.polyfit(x[-4:] - x[-1], y[-4:], 3)
        return CubicSpline(x, y, bc_type=((2, 2 * cl[1]), (2, 2 * cr[1])))

    xs = np.linspace(0, 1, 6)
    f, fp = (lambda x: np.sin(2 * x) + np.exp(x)), (lambda x: 2 * np.cos(2 * x) + np.exp(x))
    s = spline3(xs, f(xs))
    yp = s(xs, 1)
    yp0 = yp.copy()
    yp0[[0, -1]] = 0.0
    xx = np.linspace(0.01, 0.99, 99)
    sep = spline3(xs, yp)(xx)
    sep0 = spline3(xs, yp0)(xx)
    print("numpy design check, f = sin(2x) + exp(x), N=6 nodes on [0,1]:")
    print(f"  max|value' - f'|                       = {np.max(np.abs(s(xx, 1) - fp(xx))):.3e}  (interp. accuracy)")
    print(f"  max|sep.spline(yp) - value'| (no zero) = {np.max(np.abs(sep - s(xx, 1))):.3e}")
    print(f"  max|sep.spline(yp) - f'|     (no zero) = {np.max(np.abs(sep - fp(xx))):.3e}")
    print(f"  max|sep.spline(yp0) - f'| (as coded)   = {np.max(np.abs(sep0 - fp(xx))):.3e}   (|f'| ~ {np.max(fp(xx)):.2f})")


def main():
    numpy_design_check()
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(3) as ex:
        for name, st in ex.map(run, CASES):
            print(f"{name:14s} {st['status']}")
    for name in CASES:
        if name.endswith("ref"):
            continue
        ref = OUT / f"{name.split('_')[0]}_ref.npz"
        d, r = harness.load(OUT / f"{name}.npz"), harness.load(ref)
        a = d["eft_a_grid"]
        m = a > 0.01
        c = harness.compare(OUT / f"{name}.npz", ref)
        dv = np.max(np.abs(d["eft_EFTOmegaV"][m] - r["eft_EFTOmegaV"][m]))
        dp = np.max(np.abs(d["eft_EFTOmegaP"][m] - r["eft_EFTOmegaP"][m]))
        dv_lt1 = np.max(np.abs(d["eft_EFTOmegaV"][m & (a < 1)] - r["eft_EFTOmegaV"][m & (a < 1)]))
        print(f"{name:14s} vs {ref.stem}: worst Cl {c['worst_cl']:.4e}  P(k) {c['worst_pk']:.4e}  "
              f"sigma8 {c['sigma8']:.4e} | max dOmegaV {dv:.2e} (a<1: {dv_lt1:.2e})  max dOmegaP {dp:.2e} "
              f"| OmegaV(a=1) {d['eft_EFTOmegaV'][-1]:.4f} vs {r['eft_EFTOmegaV'][-1]:.4f}")
        tt = c["cl"]["total"]
        print(f"{'':14s}   total TT {tt['TT']:.3e} EE {tt['EE']:.3e} TE {tt['TE']:.3e}")
        lmax = 2500
        ct, cr = d["cl_total"][2:lmax + 1, 0], r["cl_total"][2:lmax + 1, 0]
        rel = np.abs(ct / cr - 1)
        ell = np.arange(2, lmax + 1)
        print(f"{'':14s}   TT rel diff: l<30 max {rel[ell < 30].max():.3e}, 30<=l<500 max "
              f"{rel[(ell >= 30) & (ell < 500)].max():.3e}, l>=500 max {rel[ell >= 500].max():.3e}")


if __name__ == "__main__":
    main()

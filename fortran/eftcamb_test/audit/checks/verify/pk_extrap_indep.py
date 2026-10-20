"""Independent verification of finding pk-extrap-slope-bias (camb/results.py get_matter_power_interpolator extrap_kmax).

Method (independent of the repro's code path):
  * take the raw linear/nonlinear P(k) grid from get_linear_matter_power_spectrum of a low-kmax run,
  * build the power-law extrapolation by hand with (a) the EFTCAMB secant slope mean(diff(logP[-10:]))/dlogk_last,
    (b) the upstream end slope from the last two nodes, (c) a proper secant over the last 9 intervals in log k,
  * compare with the 'truth' = log-log np.interp of the raw grid of a kmax=100/Mpc run (no CAMB interpolator),
  * check that the library interpolator of each tree (A0 = audit, B0 = upstream) reproduces the hand-built
    extrapolation on the SAME data, so any A0/B0 difference is the Python extrapolation, not the Fortran P(k),
  * test log_interp=False (non-log branch) and EFT test models read from ini (A0 only).
Run: /Users/marco/Software/miniforge3/bin/python pk_extrap_indep.py [out.json]
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else (
    harness.AUDIT_ROOT / "work" / "verify" / "pkextrap" / "pk_extrap_indep.json")
OUT.parent.mkdir(parents=True, exist_ok=True)
PARAMS = harness.PARAMS_DIR

WORKER = r"""
import sys, json, os
import numpy as np
tree = sys.argv[1]; cases = json.loads(sys.argv[2]); workdir = sys.argv[3]
sys.path.insert(0, tree)
import camb
assert os.path.realpath(camb.__file__).startswith(os.path.realpath(tree)), camb.__file__
ratios = np.array([1.5, 2.0, 5.0, 10.0])

def make_pars(case, kmax):
    if case["kind"] == "py":
        p = camb.set_params(H0=70.0, ombh2=0.0224, omch2=0.118, As=2.0e-9, ns=0.96, WantCls=False,
                            WantTransfer=True, kmax=kmax, redshifts=list(case["z"]), NonLinear=case["nl"])
        return p
    w = os.path.join(workdir, case["name"] + "_k%g.ini" % kmax)
    with open(w, "w") as f:
        f.write("DEFAULT(%s)\n" % case["ini"])
        f.write("output_root = %s\nEFTCAMB_write_background = F\n" % os.path.join(workdir, case["name"]))
        f.write("get_scalar_cls = F\nget_tensor_cls = F\nget_vector_cls = F\ndo_lensing = F\n")
        f.write("get_transfer = T\ntransfer_kmax = %g\ntransfer_num_redshifts = 1\ntransfer_redshift(1) = 0\n" % kmax)
        f.write("do_nonlinear = 0\n")
    return camb.read_ini(w, no_validate=True)

def manual(logk, logP, logk_eval):
    dl = logk[-1] - logk[-2]
    s_eft = np.mean(np.diff(logP[-10:])) / dl
    s_end = (logP[-1] - logP[-2]) / dl
    s_sec = (logP[-1] - logP[-10]) / (logk[-1] - logk[-10])
    d = logk_eval - logk[-1]
    return {"s_eft": s_eft, "s_end": s_end, "s_sec": s_sec,
            "P_eft": np.exp(logP[-1] + s_eft * d), "P_end": np.exp(logP[-1] + s_end * d),
            "P_sec": np.exp(logP[-1] + s_sec * d)}

res = {}
for case in cases:
    nlflag = case.get("nl", "NonLinear_none") != "NonLinear_none"
    out = {}
    try:
        lo = camb.get_results(make_pars(case, case.get("kmax_lo", 2.0)))
        hi = camb.get_results(make_pars(case, 100.0))
        kh, zs, P = lo.get_linear_matter_power_spectrum(nonlinear=nlflag)
        khh, zsh, Ph = hi.get_linear_matter_power_spectrum(nonlinear=nlflag)
        iz = 0  # z index of the first requested redshift (zs ordered as returned)
        zval = float(zs[iz])
        jz = int(np.argmin(np.abs(np.asarray(zsh) - zval)))
        logk = np.log(kh); logP = np.log(P[iz])
        kev = kh[-1] * ratios
        truth = np.exp(np.interp(np.log(kev), np.log(khh), np.log(Ph[jz])))
        m = manual(logk, logP, np.log(kev))
        # local end slope of the truth (hi run) at kh_max and at 10 kh_max
        lt = np.log(Ph[jz]); lkh = np.log(khh)
        def tslope(x):
            i = np.searchsorted(lkh, np.log(x)); return float((lt[i] - lt[i - 1]) / (lkh[i] - lkh[i - 1]))
        dlk = np.diff(logk[-11:])
        lib = lo.get_matter_power_interpolator(nonlinear=nlflag, extrap_kmax=30 * kh[-1], silent=True)
        P_lib = np.asarray(lib.P(zval, kev)).ravel()
        out = {"z": zval, "kh_max": float(kh[-1]), "nk": int(len(kh)), "last10_dlogk": dlk.tolist(),
               "truth_slope_at_khmax": tslope(kh[-1] * 1.0001), "truth_slope_at_10khmax": tslope(kh[-1] * 10),
               "slope_eft_formula": float(m["s_eft"]), "slope_end_upstream": float(m["s_end"]),
               "slope_true_secant": float(m["s_sec"]),
               "err_lib": (P_lib / truth - 1).tolist(),
               "err_manual_eft": (m["P_eft"] / truth - 1).tolist(),
               "err_manual_end": (m["P_end"] / truth - 1).tolist(),
               "err_manual_secant": (m["P_sec"] / truth - 1).tolist(),
               "lib_vs_manual_eft": (P_lib / m["P_eft"] - 1).tolist(),
               "lib_vs_manual_end": (P_lib / m["P_end"] - 1).tolist()}
        try:
            libn = lo.get_matter_power_interpolator(nonlinear=nlflag, extrap_kmax=30 * kh[-1], silent=True,
                                                   log_interp=False)
            Pn = np.asarray(libn.P(zval, kev)).ravel()
            out["err_lib_loginterpF"] = (Pn / truth - 1).tolist()
        except Exception as e:
            out["err_lib_loginterpF"] = "raised: " + repr(e)[:200]
        # cross spectrum delta_tot x Weyl: sign structure and extrapolation error
        try:
            kx, zx, Px = lo.get_linear_matter_power_spectrum(var1="delta_tot", var2="Weyl", nonlinear=False)
            kxh, zxh, Pxh = hi.get_linear_matter_power_spectrum(var1="delta_tot", var2="Weyl", nonlinear=False)
            out["dW_sign"] = {"all_pos": bool(np.all(Px > 0)), "all_neg": bool(np.all(Px < 0))}
            kevx = kx[-1] * ratios
            tx = -np.exp(np.interp(np.log(kevx), np.log(kxh), np.log(-Pxh[jz]))) if np.all(Pxh < 0) else None
            for li in (True, False):
                try:
                    ix = lo.get_matter_power_interpolator(var1="delta_tot", var2="Weyl", nonlinear=False,
                                                          extrap_kmax=30 * kx[-1], silent=True, log_interp=li)
                    Pxe = np.asarray(ix.P(zval, kevx)).ravel()
                    out["dW_err_loginterp_%s" % li] = (Pxe / tx - 1).tolist() if tx is not None else Pxe.tolist()
                except Exception as e:
                    out["dW_err_loginterp_%s" % li] = "raised: " + repr(e)[:200]
        except Exception as e:
            out["dW_error"] = repr(e)[:300]
    except Exception as e:
        out = {"error": repr(e)[:500]}
    res[case["name"]] = out
print("@@" + json.dumps(res))
"""


def run_tree(tree: str, cases: list[dict]) -> dict:
    workdir = OUT.parent / f"w_{tree}"
    workdir.mkdir(exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-c", WORKER, str(harness.TREES[tree]), json.dumps(cases), str(workdir)],
                           capture_output=True, text=True, env=env, cwd=workdir, timeout=1500,
                           preexec_fn=lambda: os.nice(10))
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@")]
    if not line:
        return {"_failed": p.stderr[-3000:]}
    return json.loads(line[0][2:])


def main() -> None:
    which = os.environ.get("PKX_PART", "py")
    py_cases = [
        {"kind": "py", "name": "lcdm_lin_k2", "nl": "NonLinear_none", "z": [0.0], "kmax_lo": 2.0},
        {"kind": "py", "name": "lcdm_lin_k1_z1", "nl": "NonLinear_none", "z": [1.0, 0.0], "kmax_lo": 1.0},
        {"kind": "py", "name": "lcdm_lin_k10", "nl": "NonLinear_none", "z": [0.0], "kmax_lo": 10.0},
        {"kind": "py", "name": "lcdm_nl_k2", "nl": "NonLinear_both", "z": [0.0], "kmax_lo": 2.0},
    ]
    ini_cases = [
        {"kind": "ini", "name": m, "ini": str(PARAMS / f"{m}.ini"), "kmax_lo": 2.0}
        for m in ("4_Des_F_R_LCDM_1", "2_PEFT_Omega_const_1", "5_quint_exponential", "3_RPH_KT_power_1", "5_Kmouflage_1")
    ]
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    if which == "py":
        for tree in ("A0", "B0"):
            res[tree + "_py"] = run_tree(tree, py_cases)
    else:
        res["A0_ini"] = run_tree("A0", ini_cases)
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

"""t8b check: camb/results.py get_matter_power_interpolator extrapolation (EFTCAMB MOD) vs upstream CAMB.

1) log branch: EFTCAMB replaced the upstream 2-point end slope with mean(diff(last 10)) / (last dlogk).
   Check the uniformity of the k grid over the last 10 nodes and quantify the change of P(k > kmax)
   relative to upstream for a GR model (EFTflag=0), linear and non-linear.
2) zero-crossing branch (new): check it reproduces the power-law extrapolation for a positive spectrum when
   log_interp=False, and test a spectrum with a sign change in the last 10 nodes.
Usage: python pk_extrap.py
"""

from __future__ import annotations

import sys

import numpy as np


def main(args):
    import camb

    pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2.1e-9, ns=0.965, WantCls=False,
                           WantTransfer=True, kmax=args.get("kmax", 2.0), redshifts=[2.0, 1.0, 0.5, 0.0],
                           NonLinear=args.get("nonlinear", "NonLinear_none"))
    res = camb.get_results(pars)
    nonlin = args.get("nonlinear", "NonLinear_none") != "NonLinear_none"
    khs, zs, pk = res.get_linear_matter_power_spectrum(nonlinear=nonlin)
    dlk = np.diff(np.log(khs))
    out = {"kh_max": float(khs[-1]), "last10_dlogk": dlk[-9:].tolist(),
           "dlogk_last_over_mean_last9": float(dlk[-1] / np.mean(dlk[-9:]))}
    k_test = np.array([1.0, 1.5, 2.0, 5.0, 10.0, 30.0]) * khs[-1]
    extrap = 30 * khs[-1]
    PK = res.get_matter_power_interpolator(nonlinear=nonlin, extrap_kmax=extrap, silent=True)
    P_code = PK.P(0.0, k_test)
    # upstream formula (CAMB branch `camb`): slope from the last two nodes
    logk = np.log(khs)
    iz0 = int(np.argmin(np.abs(zs)))
    out["z_order"] = list(map(float, zs))
    lp = np.log(pk[iz0])
    up_slope = (lp[-1] - lp[-2]) / (logk[-1] - logk[-2])
    eft_slope = np.mean(np.diff(lp[-10:])) / (logk[-1] - logk[-2])
    exact_slope10 = (lp[-1] - lp[-10]) / (logk[-1] - logk[-10])
    P_up = np.exp(lp[-1] + up_slope * (np.log(k_test) - logk[-1]))
    out.update({
        "z0_slope_upstream_2pt": float(up_slope),
        "z0_slope_eftcamb_mean10": float(eft_slope),
        "z0_slope_secant_10pt": float(exact_slope10),
        "k_over_kmax": (k_test / khs[-1]).tolist(),
        "P_interp_code": np.asarray(P_code).tolist(),
        "P_upstream_powerlaw": P_up.tolist(),
        "rel_code_vs_upstream": (np.asarray(P_code) / P_up - 1).tolist(),
    })
    # zero-crossing branch on a positive spectrum (log_interp=False)
    try:
        PKn = res.get_matter_power_interpolator(nonlinear=nonlin, extrap_kmax=extrap, silent=True, log_interp=False)
        out["rel_nonlog_vs_log"] = (np.asarray(PKn.P(0.0, k_test)) / np.asarray(P_code) - 1).tolist()
    except Exception as e:  # noqa: BLE001
        out["rel_nonlog_vs_log"] = []
        out["nonlog_error"] = repr(e)[:200]
    # cross spectrum that changes sign: delta_tot x v_newtonian_baryon-cdm relative velocity
    try:
        PKx = res.get_matter_power_interpolator(nonlinear=False, var1="delta_tot", var2="v_baryon_cdm",
                                                extrap_kmax=extrap, silent=True)
        kx, zx, px = res.get_linear_matter_power_spectrum(var1="delta_tot", var2="v_baryon_cdm")
        out["cross_sign_changes_last10"] = int(np.sum(np.diff(np.sign(px[iz0, -10:])) != 0))
        out["cross_islog"] = bool(PKx.islog)
        out["cross_P_beyond"] = np.asarray(PKx.P(0.0, k_test)).tolist()
        out["cross_P_last_nodes"] = px[iz0, -5:].tolist()
    except Exception as e:  # noqa: BLE001
        out["cross_error"] = repr(e)
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, WORK, run_in_tree

    allres = {}
    for tree in ("A0", "B0"):
        for nl in ("NonLinear_none", "NonLinear_both"):
            allres[f"{tree}_{nl}"] = run_in_tree(tree, HERE / "pk_extrap.py", {"nonlinear": nl, "kmax": 2.0})
    (WORK / "out" / "pk_extrap.json").write_text(json.dumps(allres, indent=1))
    for k, v in allres.items():
        if "error" in v:
            print(k, v)
            continue
        print(k, "kh_max=%.3f" % v["kh_max"], "dlogk_last/mean9=%.3f" % v["dlogk_last_over_mean_last9"],
              "slopes up/eft/secant=%.4f %.4f %.4f" % (v["z0_slope_upstream_2pt"], v["z0_slope_eftcamb_mean10"],
                                                      v["z0_slope_secant_10pt"]))
        print("   rel code vs upstream power law at k/kmax", v["k_over_kmax"], ["%.2e" % x for x in v["rel_code_vs_upstream"]])
        print("   rel nonlog vs log", ["%.2e" % x for x in v["rel_nonlog_vs_log"]])
        print("   cross:", {kk: v.get(kk) for kk in ("cross_sign_changes_last10", "cross_islog", "cross_error", "nonlog_error")})

"""t8b check: non-log (zero-crossing) extrapolation branch of get_matter_power_interpolator vs a direct high-kmax run,
for a positive spectrum (log_interp=False forced) and for a sign-changing cross spectrum (delta_cdm x v_newtonian_cdm
at low kmax is positive; we use delta_tot x Weyl which is negative definite, and a synthetic sign change).
Usage: python pk_extrap_nonlog.py TREE [TREE...]"""

import sys

import numpy as np


def main(args):
    import camb

    def pars(kmax):
        return camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2.1e-9, ns=0.965, WantCls=False,
                               WantTransfer=True, kmax=kmax, redshifts=[0.0])

    lo = camb.get_results(pars(2.0))
    hi = camb.get_results(pars(100.0))
    kh = lo.get_linear_matter_power_spectrum()[0][-1]
    k = kh * np.array([1.2, 1.5, 2.0, 5.0, 10.0])
    out = {}
    for v1, v2 in (("delta_tot", "delta_tot"), ("delta_tot", "Weyl")):
        tru = hi.get_matter_power_interpolator(nonlinear=False, var1=v1, var2=v2, silent=True, log_interp=False).P(0.0, k)
        try:
            ex = lo.get_matter_power_interpolator(nonlinear=False, var1=v1, var2=v2, extrap_kmax=30 * kh, silent=True,
                                                  log_interp=False).P(0.0, k)
            out[f"{v1}x{v2}_nonlog"] = (np.asarray(ex) / np.asarray(tru) - 1).tolist()
        except Exception as e:  # noqa: BLE001
            out[f"{v1}x{v2}_nonlog"] = repr(e)[:150]
        try:
            ex = lo.get_matter_power_interpolator(nonlinear=False, var1=v1, var2=v2, extrap_kmax=30 * kh, silent=True).P(0.0, k)
            out[f"{v1}x{v2}_log"] = (np.asarray(ex) / np.asarray(tru) - 1).tolist()
        except Exception as e:  # noqa: BLE001
            out[f"{v1}x{v2}_log"] = repr(e)[:150]
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, run_in_tree

    for t in sys.argv[1:]:
        r = run_in_tree(t, HERE / "pk_extrap_nonlog.py", {})
        print(t.split("/")[-1])
        for kk, vv in r.items():
            print("  ", kk, vv if isinstance(vv, str) else ["%+.3f" % x for x in vv])

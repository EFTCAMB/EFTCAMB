"""t8b check: behaviour of the lensing sanity check (EFTCAMB MOD in lensing.f90) for a large lensing amplitude.
Upstream raises an error ('You need to normalize realistically'); EFTCAMB silently skips lensing.
Usage: python lens_amplitude.py TREE..."""

import sys

import numpy as np


def main(args):
    import camb

    out = {}
    for As, Alens in ((2.1e-9, 1.0), (1.95e-8, 1.0), (2.1e-9, args.get("Alens_big", 300.0))):
        pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=As, ns=0.965, lmax=2500, Alens=Alens)
        try:
            res = camb.get_results(pars)
            cls = res.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)
            lens = cls["lensed_scalar"]
            unl = cls["unlensed_scalar"]
            pp = cls["lens_potential"]
            out[f"As={As},Alens={Alens}"] = {
                "error": None,
                "lensed_TT_200": float(lens[200, 0]),
                "unlensed_TT_200": float(unl[200, 0]),
                "lensed_BB_max": float(np.max(lens[2:2000, 2])),
                "lensed_all_zero": bool(np.all(lens[2:2000] == 0)),
                "PP_10_L4": float(pp[10, 0]),
            }
        except Exception as e:  # noqa: BLE001
            out[f"As={As},Alens={Alens}"] = {"error": repr(e)[:200]}
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, run_in_tree

    for t in sys.argv[1:]:
        print(t, json.dumps(run_in_tree(t, HERE / "lens_amplitude.py", {}), indent=1))

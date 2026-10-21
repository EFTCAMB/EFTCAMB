"""t8b check: for a model whose theta->H0 solve fails, scan H0 to see where calc_background_no_thermo raises
('EFTCAMB: theory unstable'), i.e. whether the failure comes from brentq trial points far from the solution."""

import sys

import numpy as np


def main(args):
    import camb

    from t8b_common import read_model

    cp = read_model(camb, args["ini"])
    out = {"H0_ini": cp.H0, "scan": {}}
    res = camb.get_background(cp, no_thermo=True)
    out["theta_at_ini"] = res.cosmomc_theta()
    for H0 in args.get("H0s", [35, 40, 41.6, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]):
        cp.H0 = H0
        try:
            camb.get_background(cp, no_thermo=True)
            out["scan"][str(H0)] = "ok"
        except Exception as e:  # noqa: BLE001
            out["scan"][str(H0)] = repr(e)[-40:]
    cp.H0 = out["H0_ini"]
    try:
        cp.set_H0_for_theta(out["theta_at_ini"], cosmomc_approx=True, theta_H0_range=(60, 80))
        out["narrow_bracket_H0"] = cp.H0
    except Exception as e:  # noqa: BLE001
        out["narrow_bracket_H0"] = repr(e)[-80:]
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, run_in_tree

    print(json.dumps(run_in_tree(sys.argv[1], HERE / "theta_unstable_scan.py", {"ini": sys.argv[2]}), indent=1))

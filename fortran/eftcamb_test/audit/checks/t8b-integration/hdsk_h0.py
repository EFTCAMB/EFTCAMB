"""t8b check: Horndeski (EFTflag=5) models: is the H0 of the EFT background equal to CP%H0 (used by CAMB for
grhocrit, Omega_x = omx h^2/h^2, halofit cosmology, theta search)? Also theta(H0) dependence."""

import sys


def main(args):
    import camb

    from t8b_common import read_model

    out = {}
    for ini in args["inis"]:
        cp = read_model(camb, ini)
        r = {"CP_H0": cp.H0}
        res = camb.get_background(cp, no_thermo=True)
        r["H_z0_dtauda"] = float(res.hubble_parameter(0.0))
        r["theta"] = res.cosmomc_theta()
        r["Omega_de_CAMB"] = res.get_Omega("de", 0.0)
        r["shooting"] = None
        for H0 in (60.0, 80.0):
            cp.H0 = H0
            try:
                res2 = camb.get_background(cp, no_thermo=True)
                r[f"H0={H0}"] = {"H_z0": float(res2.hubble_parameter(0.0)), "theta": res2.cosmomc_theta()}
            except Exception as e:  # noqa: BLE001
                r[f"H0={H0}"] = repr(e)[-60:]
        out[ini] = r
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, run_in_tree

    print(json.dumps(run_in_tree(sys.argv[1], HERE / "hdsk_h0.py", {"inis": sys.argv[2:]}, timeout=900), indent=1))

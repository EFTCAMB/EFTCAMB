"""t8b check: comoving distance (DeltaTime Romberg of dtauda) vs scipy quadrature of c/H(z) (Hofz = 1/(a^2 dtauda)).
Usage: python chi_vs_H.py TREE ini [ini...]"""

import sys

import numpy as np


def main(args):
    import camb
    from scipy.integrate import quad

    from t8b_common import read_model

    c = 299792.458
    out = {}
    for ini in args["inis"]:
        cp = read_model(camb, ini)
        res = camb.get_background(cp, no_thermo=True)
        zs = np.array([0.05, 0.2, 0.5, 1.0, 2.0, 3.0])
        chi = res.comoving_radial_distance(zs, tol=1e-10)
        chi_q = np.array([quad(lambda zz: c / res.hubble_parameter(zz), 0, z, epsabs=0, epsrel=1e-11, limit=500)[0]
                          for z in zs])
        # local derivative test at several step sizes
        d = {}
        for h in (1e-2, 1e-3, 1e-4):
            dchi = (res.comoving_radial_distance(1.0 + h, tol=1e-10) - res.comoving_radial_distance(1.0 - h, tol=1e-10)) / (2 * h)
            d[str(h)] = float(dchi * res.hubble_parameter(1.0) / c - 1)
        # smoothness of H on a fine grid around z=1
        zf = np.linspace(0.99, 1.01, 201)
        Hf = res.hubble_parameter(zf)
        d2 = np.diff(Hf, 2)
        out[ini] = {"chi_rel": (chi / chi_q - 1).tolist(), "deriv": d,
                    "H_second_diff_max_over_H": float(np.max(np.abs(d2)) / Hf[100]),
                    "H_second_diff_median_over_H": float(np.median(np.abs(d2)) / Hf[100])}
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, run_in_tree

    r = run_in_tree(sys.argv[1], HERE / "chi_vs_H.py", {"inis": sys.argv[2:]}, timeout=3000)
    print(json.dumps(r, indent=1))

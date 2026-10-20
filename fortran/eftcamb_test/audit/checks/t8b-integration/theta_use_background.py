"""t8b check: cosmomc_theta / sound horizon with EFTCAMB_use_background=T vs F (Horndeski, EFTflag=5).

CAMB defines cosmomc_theta with the approximate sound speed R = 3e4 a ombh2 (dsound_da_approx), while the
EFT tabulated r_sound uses the exact R = 3 rho_b / (4 rho_g). Recompute both by direct quadrature of the
same H(z) to see which definition each code path returns.
Usage: python theta_use_background.py [tree] [ini]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main(args):
    import camb
    from scipy.integrate import quad

    from t8b_common import PARAMS, WORK

    out = {}
    for ub in ("F", "T"):
        ini = WORK / "out" / f"theta_ub{ub}.ini"
        ini.write_text(f"DEFAULT({PARAMS / args['ini']})\nEFTCAMB_use_background = {ub}\nfeedback_level = 0\n")
        cp = camb.read_ini(str(ini), no_validate=True)
        res = camb.get_background(cp)
        th_mc = res.cosmomc_theta()
        der = res.get_derived_params()
        zstar_fit = None
        # CAMB's zstar fitting formula used in CosmomcTheta
        ombh2 = cp.ombh2
        omdmh2 = cp.omch2 + cp.omnuh2
        zstar_fit = 1048 * (1 + 0.00124 * ombh2**-0.738) * (
            1 + (0.0783 * ombh2**-0.238 / (1 + 39.5 * ombh2**0.763)) * (omdmh2 + ombh2) ** (0.560 / (1 + 21.1 * ombh2**1.81))
        )
        astar = 1 / (1 + zstar_fit)
        grhog_over_grhob = None

        def dtauda(a):
            return 1.0 / (a * a * float(res.h_of_z(1 / a - 1)))

        def cs_approx(a):
            return 1 / np.sqrt(3 * (1 + 3.0e4 * a * ombh2))

        # exact: R = 3 rho_b/(4 rho_g); omega_gamma h^2 from TCMB
        omgh2 = 2.4728e-5 * (cp.TCMB / 2.7255) ** 4 * (2.7255 / 2.7255) ** 0
        dens = res.get_background_densities([1.0], vars=["photon", "baryon"])
        Rfac = 3 * dens["baryon"][0] / (4 * dens["photon"][0])  # R = Rfac * a

        def cs_exact(a):
            return 1 / np.sqrt(3 * (1 + Rfac * a))

        lo = 1e-8
        pts = list(np.logspace(-8, np.log10(astar), 40))
        rs_approx = sum(quad(lambda a: dtauda(a) * cs_approx(a), pts[i], pts[i + 1], epsrel=1e-10)[0] for i in range(39))
        rs_exact = sum(quad(lambda a: dtauda(a) * cs_exact(a), pts[i], pts[i + 1], epsrel=1e-10)[0] for i in range(39))
        DA = float(res.angular_diameter_distance(zstar_fit)) / astar
        out[ub] = {
            "H0": cp.H0,
            "cosmomc_theta_code": th_mc,
            "theta_approx_quadrature": rs_approx / DA,
            "theta_exact_quadrature": rs_exact / DA,
            "rstar_derived": der["rstar"],
            "thetastar_derived": der["thetastar"],
            "zstar": der["zstar"],
            "age": der["age"],
            "rdrag": der["rdrag"],
            "H(z=1)": float(res.hubble_parameter(1.0)),
            "tau0": float(res.conformal_time(0.0)),
        }
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, run_in_tree

    tree = sys.argv[1] if len(sys.argv) > 1 else "A0"
    ini = sys.argv[2] if len(sys.argv) > 2 else "5_hdsk_TG_1.ini"
    print(json.dumps(run_in_tree(tree, HERE / "theta_use_background.py", {"ini": ini}), indent=1))

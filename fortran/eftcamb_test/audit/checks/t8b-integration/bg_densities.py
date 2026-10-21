"""t8b check: are CAMB v2 background-density APIs routed through the EFT background?

Compares H(z) from results.hubble_parameter (dtauda -> EFT) with H(z) reconstructed from
results.get_background_densities (Fortran CAMBdata_GetBackgroundDensities, uses CP%DarkEnergy),
and Omega_de from get_Omega with 1 - (non-DE densities)/H^2.
Usage: python bg_densities.py [tree] [ini ...]
"""

from __future__ import annotations

import sys

import numpy as np

C_KMS = 299792.458


def main(args):
    import camb
    from t8b_common import read_model

    out = {}
    for ini in args["inis"]:
        cp = read_model(camb, ini)
        res = camb.get_background(cp)
        z = np.array([0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0, 1000.0])
        a = 1 / (1 + z)
        H = np.asarray(res.hubble_parameter(z))
        d = res.get_background_densities(a, vars=["tot", "de", "cdm", "baryon", "photon", "neutrino", "nu", "K"])
        H_dens = np.sqrt(d["tot"] / 3) / a**2 * C_KMS
        nonde = d["tot"] - d["de"]
        om_de_true = 1 - nonde / (3 * (H / C_KMS) ** 2 * a**4)
        om_de_api = np.asarray(res.get_Omega("de", z))
        rho_de, w_de = res.get_dark_energy_rho_w(a)
        out[ini] = {
            "EFTflag": int(cp.EFTCAMB.EFTflag) if hasattr(cp, "EFTCAMB") else None,
            "z": z.tolist(),
            "H_dtauda": H.tolist(),
            "H_from_densities": H_dens.tolist(),
            "max_rel_H_mismatch": float(np.max(np.abs(H_dens / H - 1))),
            "Omega_de_api": om_de_api.tolist(),
            "Omega_de_from_H": om_de_true.tolist(),
            "max_abs_Omega_de_mismatch": float(np.max(np.abs(om_de_api - om_de_true))),
            "w_de_api": np.asarray(w_de).tolist(),
        }
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, run_in_tree

    tree = sys.argv[1] if len(sys.argv) > 1 else "A0"
    inis = sys.argv[2:] or ["1_EFT_GR.ini", "2_PEFT_wDE_wCDM_1.ini", "4_Des_mc5e_CPL_1.ini", "5_Kmouflage_1.ini"]
    print(json.dumps(run_in_tree(tree, HERE / "bg_densities.py", {"inis": inis}), indent=1))

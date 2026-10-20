"""t8b check: are CAMB v2 background outputs mutually consistent with the EFT expansion history (H(z) from dtauda)?
For each model: (a) chi(z) vs quadrature of c/H(z); (b) age vs quadrature of 1/((1+z)H); (c) r_s(z_drag) vs quadrature of
c_s/H with the exact baryon loading; (d) BAO rs/D_V vs recomputation; (e) tau0 vs conformal_time(0).
Usage: python bg_routing.py TREE"""

import sys

import numpy as np


def main(args):
    import camb
    from scipy.integrate import quad

    from t8b_common import read_model

    c = 299792.458
    out = {}
    for ini, ub in args["cases"]:
        cp = read_model(camb, ini)
        if ub is not None:
            cp.EFTCAMB.EFTCAMB_use_background = ub
        res = camb.get_background(cp)
        z = np.array([0.1, 0.5, 1.0, 2.0, 5.0])
        chi = res.comoving_radial_distance(z, tol=1e-10)
        chi_q = np.array([quad(lambda zz: c / res.hubble_parameter(zz), 0, zz1, epsabs=0, epsrel=1e-11, limit=500)[0]
                          for zz1 in z])
        a_err = float(np.max(np.abs(chi / chi_q - 1)))
        Mpc_km = 3.0856775814913673e19
        Gyr_s = 3.15576e16
        age_q = quad(lambda zz: 1 / ((1 + zz) * res.hubble_parameter(zz)), 0, np.inf, limit=400)[0] * Mpc_km / Gyr_s
        age = res.physical_time(0.0)
        der = res.get_derived_params()
        zd = der["zdrag"]
        Tcmb = cp.TCMB
        ombh2 = cp.ombh2
        # photon density omega_g h^2
        ogh2 = 4.48162687719e-7 * Tcmb**4

        def cs(zz):
            R = 0.75 * ombh2 / ogh2 / (1 + zz)
            return c / np.sqrt(3 * (1 + R))

        rs_q = quad(lambda zz: cs(zz) / res.hubble_parameter(zz), zd, np.inf, limit=400, epsrel=1e-9)[0]
        rs = res.sound_horizon(zd)
        zb = np.array([0.38, 0.61])
        DV = (res.angular_diameter_distance(zb) ** 2 * (1 + zb) ** 2 * c * zb / res.hubble_parameter(zb)) ** (1 / 3)
        bao = res.get_BAO(zb, cp)[:, 0]
        out[f"{ini}|ub={ub}"] = {
            "chi_vs_quad_c_over_H": a_err,
            "age_rel": age / age_q - 1,
            "rdrag_code_vs_quad": rs / rs_q - 1,
            "rdrag_derived_vs_sound_horizon": der["rdrag"] / rs - 1,
            "bao_rs_DV_rel": float(np.max(np.abs(bao / (der["rdrag"] / DV) - 1))),
            "tau0_rel": res.conformal_time(0.0) / res.tau0 - 1,
        }
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, WORK, run_in_tree

    cases = [["1_EFT_GR.ini", None], ["2_PEFT_wDE_wCDM_1.ini", None], ["4_Des_mc5e_CPL_1.ini", None],
             ["5_Kmouflage_1.ini", None], ["5_quint_exponential.ini", None], ["3_RPH_MK_OmegaDE_1.ini", None],
             ["5_hdsk_TG_1.ini", False], ["5_hdsk_TG_1.ini", True]]
    r = run_in_tree(sys.argv[1], HERE / "bg_routing.py", {"cases": cases}, timeout=3000)
    (WORK / "out" / "bg_routing.json").write_text(json.dumps(r, indent=1))
    if "error" in r and "stdout_tail" in r:
        print(r)
    for k, v in r.items():
        print(f"{k:34s} " + " ".join(f"{kk}={vv:+.1e}" for kk, vv in v.items()))

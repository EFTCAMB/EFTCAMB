"""Check 8: background getters with EFT on: get_background_densities, get_dark_energy_rho_w, get_Omega, get_age,
physical_time, conformal_time, sound horizon, get_zre_from_tau, get_background_redshift_evolution.

Internal consistency with H(z) (hubble_parameter uses the EFT-aware dtauda):
  Friedmann: sum_i 8 pi G a^4 rho_i (get_background_densities 'tot') == 3 a^4 (H/c)^2
  age: get_age == physical_time(0) == int_0^1 da / (a H(a))   (numerical, from hubble_parameter)
  tau0 == int_0^1 da / (a^2 H)
  get_zre_from_tau(params, tau_ini) == Reion.redshift of the run that used tau_ini
"""

from __future__ import annotations

import json
import sys

import numpy as np
from scipy import integrate
from t8c_lib import WORK, load, model_ini, run_many

OUT = WORK / "c8"
C_KMS = 299792.458
OV = {"get_tensor_cls": "F", "get_scalar_cls": "F", "get_transfer": "F", "feedback_level": "0",
      "num_redshiftwindows": "0"}
CASES = {
    "GR": ("1_EFT_GR", {}),
    "GRw08": ("1_EFT_GR", {"w": "-0.8"}),
    "PEFTw08": ("2_PEFT_wDE_wCDM_1", {}),
    "fRw09": ("4_Des_F_R_wCDM_1", {}),
    "quint": ("5_quint_exponential", {}),
    "Kmouf": ("5_Kmouflage_1", {}),
    "hdsk": ("5_hdsk_TG_1", {}),
    "hdsk_usebg": ("5_hdsk_TG_1", {"EFTCAMB_use_background": "T"}),
}


def jobs():
    return [dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{tag}.npz", overrides=OV | ov,
                 args={"outputs": ["background_only", "background", "zre"], "zre_taus": [0.09, 0.054]})
            for tag, (m, ov) in CASES.items()]


def analyse(d):
    z = d["bg_z"]
    a = 1 / (1 + z)
    H = d["bg_H"]
    tot = d["bg_dens"][:, 0]
    fried = tot / (3 * a**4 * (H / C_KMS) ** 2) - 1
    m = d["meta"]
    # numerical age and tau0 from H(z) on a fine grid (hubble_parameter only, via a spline in ln a)
    res = {"friedmann_maxrel": float(np.max(np.abs(fried))), "friedmann_z0": float(fried[0]),
           "de_w_z0": float(d["bg_de_w"][0]), "Omega_de_z0": float(d["bg_Omega_de"][0]),
           "age_get_age": m["age_get_age"], "age_physical_time": m["age_results"],
           "age_derived": m.get("derived", {}).get("age"), "tau0": m["tau0"],
           "zre_from_tau_0.09": m.get("zre_from_tau_0.09"), "reion_redshift": m.get("reion_redshift"),
           "zrei_derived": m.get("reion_tau_derived_zre"), "errors": m["errors"]}
    if "eft_grhov_t" in d:
        res["eft_fields"] = True
    return res


def numerical(tag):
    """Age and tau0 from hubble_parameter on a dense grid (needs a separate run with dense z: use bg_H spline)."""
    d = load(OUT / f"{tag}.npz")
    z = d["bg_z"]
    H = d["bg_H"]
    lna = -np.log1p(z)
    o = np.argsort(lna)
    lna, H = lna[o], H[o]
    # integrand for age: dt = da/(aH) = dlna / H ; for tau: dtau = dlna /(a H)
    from scipy.interpolate import CubicSpline

    s = CubicSpline(lna, np.log(H))
    x = np.linspace(lna[0], 0, 20001)
    Hx = np.exp(s(x))
    age_tail = integrate.simpson(1 / Hx, x=x)  # from z=1e4 to 0 in Mpc s/km
    tau_tail = integrate.simpson(1 / (np.exp(x) * Hx), x=x) * C_KMS
    Mpc_km = 3.0856775814913673e19
    Gyr = 3600 * 24 * 365.25e9
    return {"age_from_H_z<1e4_Gyr": float(age_tail * Mpc_km / Gyr), "tau_from_H_z<1e4": float(tau_tail)}


def report():
    res = {}
    for tag in CASES:
        try:
            d = load(OUT / f"{tag}.npz")
        except FileNotFoundError:
            st = json.loads((OUT / f"{tag}.json").read_text())
            res[tag] = st["status"] + " " + str(st.get("log_tail", [])[-2:])
            continue
        res[tag] = analyse(d) | numerical(tag)
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        for r in run_many(jobs(), nproc=int(sys.argv[2]) if len(sys.argv) > 2 else 4):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report()

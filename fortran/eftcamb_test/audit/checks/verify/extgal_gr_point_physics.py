"""Verifier (physics lens) for t3-extgal-gr-point-background.

Runs in-process on the A0 tree (one CPU slot). Independent expectation: the Extended Galileon tracker equation
    E^(2+S) = Omega_phi0 + Omega_m(a) E^S,   Omega_m(a) = non-DE density / (3 H0^2)   (incl. massive nu)
with Omega_m(a) and Omega_phi0 taken from a plain LCDM run (EFTflag=0) at the same cosmology.
"""

import os
import sys
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

sys.path.insert(0, str(harness.AUDIT_ROOT / "trees" / "A0"))
import camb  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

INI = Path(__file__).with_name("extgal_fixedpath.ini")
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
FAST = dict(get_tensor_cls=False, do_lensing=False)
C = 299792.458


WORK = harness.AUDIT_ROOT / "work" / "verify" / "extgal_gr_point"


def run(ini, **ov):
    WORK.mkdir(parents=True, exist_ok=True)
    name = Path(ini).stem + "_" + "_".join(f"{k}{v}" for k, v in ov.items()) + ".ini"
    w = WORK / name
    w.write_text("".join(f"{k} = {v}\n" for k, v in ov.items()) + f"DEFAULT({Path(ini).resolve()})\n")
    cp = camb.read_ini(str(w), no_validate=True)
    cp.WantTensors = False
    cp.InitPower.r = 0.0
    cp.max_l = 1500
    cp.WantTransfer = False
    return cp, camb.get_results(cp)


def tracker_E(a, S, omphi0, omtot):
    if S == 0:
        return np.sqrt(omphi0 + omtot)
    f = lambda E: omphi0 + omtot * E**S - E ** (2 + S)  # noqa: E731
    return brentq(f, 1e-6, 1e9, xtol=1e-15, rtol=1e-14)


def main():
    zs = np.array([0.0, 1.0, 3.0, 3.9, 4.1, 5.0, 9.0, 30.0, 100.0])
    a_list = [0.1, 0.19, 0.21, 0.5, 1.0]
    with harness.slot(1):
        for mnu_label, mnu in [x for x in (("mnu_base", None), ("mnu0", 0.0)) if x[0] in os.environ.get("MNU", "mnu_base,mnu0")]:
            ov = {} if mnu is None else {"omnuh2": mnu}
            cpL, rL = run(GR, **ov)
            H0 = cpL.H0
            dens = rL.get_background_densities(np.array([1.0 / (1 + z) for z in zs] + a_list))
            tot, de = dens["tot"], dens["de"]
            omphi0 = rL.get_background_densities(1.0)["de"][0] / rL.get_background_densities(1.0)["tot"][0]
            HL = rL.hubble_parameter(zs)
            om_nonde_z = (HL / H0) ** 2 * (1 - de[: len(zs)] / tot[: len(zs)])
            HL_a = rL.hubble_parameter(1 / np.array(a_list) - 1)
            om_nonde_a = (HL_a / H0) ** 2 * (1 - de[len(zs):] / tot[len(zs):])
            print(f"=== {mnu_label}: omnuh2={cpL.omnuh2:.4f}, Omega_phi0(LCDM)={omphi0:.6f}")
            for S, q, skip in [c for i, c in enumerate(((0.0, 100.0, True), (2.0, 0.35, False), (2.0, 0.35, True))) if str(i) in os.environ.get("CASES", "012")]:
                try:
                    cpE, rE = run(INI, **ov, ExtendedGalileon_S=S, ExtendedGalileon_q=q,
                                  EFTCAMB_skip_stability="T" if skip else "F")
                except Exception as e:  # noqa: BLE001
                    print(f"  S={S} q={q} skipstab={skip}: get_results FAILED: {e!r}")
                    continue
                HE = rE.hubble_parameter(zs)
                Eexp = np.array([tracker_E(1 / (1 + z), S, omphi0, o) for z, o in zip(zs, om_nonde_z)])
                print(f"  S={S} q={q} skipstab={skip}: CAMB H(z)/H_expected - 1:")
                print("     " + "  ".join(f"z={z:g}:{h / (H0 * e) - 1:+.2e}" for z, h, e in zip(zs, HE, Eexp)))
                _, arr = cpE.EFTCAMB.get_eft_functions(rE, np.array(a_list))
                ad = arr["adotoa"]
                Eexp_a = np.array([tracker_E(a, S, omphi0, o) for a, o in zip(a_list, om_nonde_a)])
                print("     get_eft_functions adotoa/(a H_exp) - 1: " + "  ".join(
                    f"a={a:g}:{x / (a * H0 * e / C) - 1:+.2e}" for a, x, e in zip(a_list, ad, Eexp_a)))
                cls = rE.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"]
                print(f"     C_l any NaN: {bool(np.isnan(cls).any())};  theta*={rE.get_derived_params()['thetastar']:.6f}"
                      f" (LCDM {rL.get_derived_params()['thetastar']:.6f})")


if __name__ == "__main__":
    main()

"""Repro (t3-gr): Extended Galileon (EFTflag=4, FullMappingEFTmodel=7): wrong expansion histories (CAMB background
without dark energy at z > 4, get_eft_functions adotoa without massive neutrinos at a >= 0.2) and NaN spectra at the
LCDM point.

At S = 0 (B = S q = 0) the model is exactly LCDM: EFTc = Gamma_i = 0, EFTLambda = -3 a^2 H0^2 Omega_phi0, and the
tracker equation E^(2+S) = Omega_phi0 + Omega_tot E^S gives E^2 = Omega_phi0 + Omega_tot.
EFTCAMBExtendedGalileonComputeAdotoa (08f_full_models/008p7_Extended_Galileon.f90:435-497):
  * a < 0.2: adotoa = sqrt((grhom_t + 2c - Lambda)/3/(...)) from the cache. Inside CAMB's dtauda
    (EFTCAMBFullModelComputeDtauda path) only grhoa2 is filled (grhom_t = c = Lambda = 0), so adotoa = 0, the
    background EFT functions return early (a*adotoa == 0), c = Lambda = 0, and CAMB's H(z) (distances, thermal
    history, results.hubble_parameter) has no dark energy for z > 4: H/H_LCDM - 1 = -0.9% at z = 4.1. This holds
    for every S, not only at the LCDM point.
  * a >= 0.2: the tracker equation uses Omega_tot with eft_cache%grhonu_tot, which is filled only AFTER
    compute_adotoa (equations.f90 derivs: lagged by one call; 11_EFTCAMB_python.f90 EFTCAMB_GetEFTFunctions: still 0),
    so massive neutrinos are dropped: get_eft_functions adotoa is low by ~Omega_nu/2 (base_params omnuh2 = 0.01:
    Omega_nu = 0.0204, -1.0% at a = 1, -3.4% at a = 0.21); with omnuh2 = 0 the error vanishes (control run).
    (The bisection stop test `sqrt(fsoluction**2) > 10**(-1)` is integer 10**(-1) = 0, so it always runs 50 halvings;
    harmless, but not what was intended.)
At S = 0 the stability check finds NaNs in the cache at a = 1e-10 ("theory unstable"); with the check skipped the
returned C_l are NaN.
Expected (S=0): H(z) and adotoa equal LCDM to ~1e-6 and finite C_l equal to EFTflag=0.
Run: /Users/marco/Software/miniforge3/bin/python t3-extgal-gr-point-background.py   (5 fast CAMB runs, <1 min)
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

INI = REPO / "fortran" / "eftcamb_test" / "audit" / "inis" / "008p7-ext-galileon.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
OUT = harness.AUDIT_ROOT / "work" / "t3-gr" / "repro_extgal"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "200", "get_transfer": "F",
        "num_redshiftwindows": "0"}
C = 299792.458


def job(args):
    name, ini, ov = args
    out = OUT / f"{name}.npz"
    st = harness.run_model("A0", ini, out, overrides={**FAST, **ov}, threads=1, timeout=600)
    return name, st, out


def main() -> None:
    jobs = [("lcdm", GR, {}),
            ("extgal_S0", INI, {"ExtendedGalileon_S": "0", "ExtendedGalileon_q": "100"}),
            ("extgal_S0_skipstab", INI, {"ExtendedGalileon_S": "0", "ExtendedGalileon_q": "100",
                                         "EFTCAMB_skip_stability": "T"}),
            ("lcdm_nu0", GR, {"omnuh2": "0.0"}),
            ("extgal_S0_skipstab_nu0", INI, {"ExtendedGalileon_S": "0", "ExtendedGalileon_q": "100",
                                             "EFTCAMB_skip_stability": "T", "omnuh2": "0.0"})]
    with ThreadPoolExecutor(4) as ex:
        res = {n: (st, out) for n, st, out in ex.map(job, jobs)}
    st, _ = res["extgal_S0"]
    print(f"S=0 with stability check: status={st['status']}",
          [ln.strip() for ln in st.get("log_tail", []) if "unstable" in ln or "Nan in" in ln][-2:])
    st, out = res["extgal_S0_skipstab"]
    if st["status"] != "ok":
        print("S=0 with stability skipped:", st["status"])
        return
    d, r = harness.load(out), harness.load(res["lcdm"][1])
    print("S=0, stability skipped: any NaN in C_l:", bool(np.isnan(d["cl_total"]).any()))
    print("CAMB background H(z)/H_LCDM(z) - 1 (dtauda):")
    for z in (1.0, 3.0, 3.9, 4.1, 5.0, 9.0, 30.0):
        print(f"   z={z:5.1f}: {np.interp(z, d['bg_z'], d['bg_H']) / np.interp(z, r['bg_z'], r['bg_H']) - 1:+.2e}")
    print("get_eft_functions adotoa / (a H_LCDM / c) - 1   [omnuh2 = 0.01 | control omnuh2 = 0]:")
    a = d["eft_a_grid"]
    d0, r0 = harness.load(res["extgal_S0_skipstab_nu0"][1]), harness.load(res["lcdm_nu0"][1])
    for av in (0.1, 0.19, 0.21, 0.5, 1.0):
        i = int(np.argmin(abs(a - av)))
        hl = np.interp(1 / a[i] - 1, r["bg_z"], r["bg_H"])
        hl0 = np.interp(1 / a[i] - 1, r0["bg_z"], r0["bg_H"])
        print(f"   a={a[i]:.3f}: {d['eft_adotoa'][i] / (a[i] * hl / C) - 1:+.2e} | "
              f"{d0['eft_adotoa'][i] / (a[i] * hl0 / C) - 1:+.2e}")
    print("expected: status ok, finite C_l, all ratios ~0 (<1e-6)")


if __name__ == "__main__":
    main()

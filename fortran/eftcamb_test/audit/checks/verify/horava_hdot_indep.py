"""Independent check of finding horava-hdot-inconsistent (verify track).

Different code path from the repro: the true dH/dtau is taken from CAMB's *background* (results.h_of_z, which goes
through the model's compute_dtauda), and the Hdot/Hdotdot under test are the values held in the EFT timestep cache
*during the perturbation evolution* (EFTCAMB_quantities_evolution -> equations.f90 derivs), not the python
get_eft_functions bridge.  F is computed from the ini parameters, not fitted.

Usage: python horava_hdot_indep.py [tree] [out_json]
"""

import os

os.environ["OMP_NUM_THREADS"] = "1"  # before numpy/camb load libgomp
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else AUDIT / "trees" / "A0"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else AUDIT / "work" / "verify" / f"horava_hdot_indep_{TREE.name}.json"
REPO = Path(__file__).resolve().parents[5]
PARAMS = REPO / "fortran" / "eftcamb_test" / "parameters"
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, str(TREE))
sys.path.insert(1, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
import camb  # noqa: E402

assert Path(camb.__file__).resolve().is_relative_to(TREE.resolve()), camb.__file__

# 8th-order central first-derivative stencil
W8 = np.array([1 / 280, -4 / 105, 1 / 5, -4 / 5, 0.0, 4 / 5, -1 / 5, 4 / 105, -1 / 280])
OFFS = np.arange(-4, 5)


def F_of(p: dict) -> float:
    eta = float(p.get("Horava_eta", 0))
    lam = float(p.get("Horava_lambda", 0))
    xi = 0.5 * eta if p.get("HoravaSolarSystem", "F").strip().upper().startswith("T") else float(p.get("Horava_xi", 0))
    return (2 + 2 * xi - eta) / (3 * lam + 2)


def adotoa_bg(res, a):
    a = np.asarray(a, dtype=float)
    return a * res.h_of_z(1 / a - 1)  # conformal Hubble rate [1/Mpc] from CAMB background (dtauda)


def run(model: str, a_c: np.ndarray, h: float = 0.01, k: float = 0.01) -> dict:
    tmp = Path(tempfile.mkdtemp(prefix="vf_horava_", dir=AUDIT / "work" / "verify"))
    ini = tmp / "w.ini"
    ini.write_text(f"DEFAULT({PARAMS / model}.ini)\noutput_root = {tmp}/out\nEFTCAMB_write_background = F\n"
                   "EFTCAMB_feedback_level = 0\nfeedback_level = 0\nget_transfer = F\n")
    p = harness.parse_ini(PARAMS / f"{model}.ini")
    F = F_of(p)
    cp = camb.read_ini(str(ini), no_validate=True)
    res = camb.get_results(cp)
    # --- true derivatives from the background, nested 8th-order FD in ln a
    H = adotoa_bg(res, a_c)
    def hdot_true(ac):
        aa = ac[:, None] * np.exp(h * OFFS[None, :])
        Hs = adotoa_bg(res, aa.ravel()).reshape(aa.shape)
        return Hs[:, 4] * (Hs @ W8) / h
    Hdot_true = hdot_true(a_c)
    h2 = 0.02
    aa2 = a_c[:, None] * np.exp(h2 * OFFS[None, :])
    Hd_s = hdot_true(aa2.ravel()).reshape(aa2.shape)
    Hdotdot_true = H * (Hd_s @ W8) / h2
    # --- cache values held during the perturbation evolution (equations.f90 path)
    fields, arr = cp.EFTCAMB.get_scale_evolution(res, [k], a_c)
    # arr shape (nk, ntimes)
    ev = {f: np.asarray(arr[f])[0] for f in arr.dtype.names}
    Hdot_cons = -0.5 * (ev["adotoa"] ** 2 + F * (ev["gpresm_t"] + ev["gpiv_t"]))
    Hdot_code_formula = np.sqrt(F) * (-0.5 * (ev["adotoa"] ** 2 + ev["gpresm_t"] + ev["gpiv_t"]))
    out = {"model": model, "F": F, "a": a_c.tolist(), "k": k,
           "adotoa_bg": H.tolist(), "adotoa_cache": ev["adotoa"].tolist(),
           "rel_adotoa_cache_vs_bg": (np.abs(ev["adotoa"] - H) / np.abs(H)).tolist(),
           "Hdot_true_bgFD": Hdot_true.tolist(), "Hdot_cache": ev["Hdot"].tolist(),
           "Hdot_consistent_formula": Hdot_cons.tolist(), "Hdot_code_formula_recomputed": Hdot_code_formula.tolist(),
           "rel_Hdot_cache_vs_true": (np.abs(ev["Hdot"] - Hdot_true) / np.abs(Hdot_true)).tolist(),
           "rel_Hdot_consistent_vs_true": (np.abs(Hdot_cons - Hdot_true) / np.abs(Hdot_true)).tolist(),
           "rel_Hdot_cache_vs_codeformula": (np.abs(ev["Hdot"] - Hdot_code_formula) / np.abs(Hdot_code_formula)).tolist(),
           "Hdotdot_true_bgFD": Hdotdot_true.tolist(), "Hdotdot_cache": ev["Hdotdot"].tolist(),
           "rel_Hdotdot_cache_vs_true": (np.abs(ev["Hdotdot"] - Hdotdot_true) / np.abs(Hdotdot_true)).tolist(),
           "Hdotdotdot_cache": ev["Hdotdotdot"].tolist(),
           "EFTLambda_cache": ev["EFTLambda"].tolist(), "EFTc_cache": ev["EFTc"].tolist()}
    # EFT functions recomputed with the true Hdot (mapping formulas of 008p1, lines 335-349)
    eta = float(p.get("Horava_eta", 0)); lam = float(p.get("Horava_lambda", 0))
    xi = 0.5 * eta if p.get("HoravaSolarSystem", "F").strip().upper().startswith("T") else float(p.get("Horava_xi", 0))
    D = 2 + 2 * xi - eta
    c_true = -(2 * xi - 3 * lam) * (Hdot_true - H**2) / D
    L_true = -ev["grhov_t"] + 2 * (3 * lam - 2 * xi) * (0.5 * H**2 + Hdot_true) / D
    out["rel_EFTc_cache_vs_trueHdot"] = (np.abs(ev["EFTc"] - c_true) / np.abs(c_true)).tolist()
    out["rel_EFTLambda_cache_vs_trueHdot"] = (np.abs(ev["EFTLambda"] - L_true) / np.abs(L_true)).tolist()
    return out


def main():
    a_c = np.array([0.02, 0.1, 0.3, 0.5, 0.9])
    allout = []
    with harness.slot(1):
        for model in ("5_Horava_eta_4", "5_Horava_SolSyst_3", "5_Horava_eta_1"):
            o = run(model, a_c)
            allout.append(o)
            print(f"== {model}  F={o['F']:.6f}  tree={TREE.name}")
            print("     a     adotoa cache/bg-1   Hdot_cache      Hdot_true(bgFD)  rel(cache)  rel(consistent)  "
                  "rel(Hdd cache)  Hddd_cache  rel EFTc  rel EFTLambda")
            for i, a in enumerate(a_c):
                print(f"  {a:5.2f}   {o['rel_adotoa_cache_vs_bg'][i]:.2e}     {o['Hdot_cache'][i]: .6e}  "
                      f"{o['Hdot_true_bgFD'][i]: .6e}   {o['rel_Hdot_cache_vs_true'][i]:.3e}   "
                      f"{o['rel_Hdot_consistent_vs_true'][i]:.3e}       {o['rel_Hdotdot_cache_vs_true'][i]:.3e}     "
                      f"{o['Hdotdotdot_cache'][i]: .2e}  {o['rel_EFTc_cache_vs_trueHdot'][i]:.2e}  "
                      f"{o['rel_EFTLambda_cache_vs_trueHdot'][i]:.2e}")
            print(f"  cache Hdot == sqrt(F)*GR-formula: max rel {max(o['rel_Hdot_cache_vs_codeformula']):.2e}")
    OUT.write_text(json.dumps(allout, indent=1))
    print("wrote", OUT)


if __name__ == "__main__":
    main()

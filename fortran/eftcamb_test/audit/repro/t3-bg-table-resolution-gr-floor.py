"""Repro (t3-gr): K-mouflage / K-mimic (008p3), ADE (008p2) and quintessence (008p4) do not converge to LCDM as
their coupling -> 0: the default background tables are too coarse.

Both models tabulate EFTOmega, EFTLambda, EFTc, EFTgamma1 on model_background_num_points points equispaced in ln a
(K-mouflage: 1000 points from a=1e-9, ADE: 1500 from a=1e-10) and read them with linear interpolation
(equispaced_linear_interpolation_1D). With d ln a = 0.021 (0.015) the interpolation error of EFTLambda ~ a^2 is
~(d ln a)^2/8 * 4 ~ 2e-4 (1e-4), i.e. an H(z) error ~Omega_DE * 1e-4 at z < 1, and C_l errors above 1e-4, even at
eps2_0 = -1e-8 / f_ac = 1e-8 where the models are LCDM to 1e-8. Ten times more points removes the floor (error falls
~100x, as expected for linear interpolation). Quintessence (1100 points from a=1e-8 to 1.1, exponential potential with
EFT_lambda = 0.003, physical deviation ~0.15 lambda^2 ~ 1e-6) shows the same ~1e-4 floor.
Expected: deviation from EFTflag=0 -> ~1e-6 as the coupling -> 0. Observed: printed.
Run: /Users/marco/Software/miniforge3/bin/python t3-bg-table-resolution-gr-floor.py   (9 CAMB runs, ~2-3 min)
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

P = harness.PARAMS_DIR
OUT = harness.AUDIT_ROOT / "work" / "t3-gr" / "repro_bgtable"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "2000", "num_redshiftwindows": "0"}
# ADE (turn-on a = 1e-10) needs the test-ini CMB settings: reduced settings trip the 'EFTCAMB starts before thermo
# tauminn' guard (see t3-eft-turnon-before-tauminn-lmax.py)
ADE_SET = {"get_tensor_cls": "F"}
JOBS = {
    "GR_nu0": (P / "1_EFT_GR.ini", {"omnuh2": "0.0"}, {"Accuracy.IntTolBoost": 2.0}),
    "GR": (P / "1_EFT_GR.ini", {"ADE": "1"}, {"Accuracy.IntTolBoost": 2.0}),
    "GR_fast": (P / "1_EFT_GR.ini", {}, {"Accuracy.IntTolBoost": 2.0}),
    "Kmouflage_eps2_1e-8": (P / "5_Kmouflage_1.ini", {"eps2_0": "-1e-08"}, {}),
    "Kmouflage_eps2_1e-8_10xpoints": (P / "5_Kmouflage_1.ini", {"eps2_0": "-1e-08",
                                                                "model_background_num_points": "10000"}, {}),
    "ADE_fac_1e-8": (P / "5_ADE_1.ini", {"ADE": "1", "f_ac": "1e-08"}, {}),
    "ADE_fac_1e-8_10xpoints": (P / "5_ADE_1.ini", {"ADE": "1", "f_ac": "1e-08", "model_background_num_points": "15000"},
                               {}),
    "quint_exp_lambda_3e-3": (P / "5_quint_exponential.ini", {"EFT_lambda": "0.003"}, {}),
    "quint_exp_lambda_3e-3_10xpoints": (P / "5_quint_exponential.ini", {"EFT_lambda": "0.003",
                                                                        "model_background_num_points": "11000"}, {}),
}


def job(item):
    name, (ini, ov, ps) = item
    out = OUT / f"{name}.npz"
    base = ADE_SET if ov.pop("ADE", None) else FAST
    st = harness.run_model("A0", ini, out, overrides={**base, **ov}, python_set=ps, threads=1, timeout=900)
    return name, st, out


def main() -> None:
    with ThreadPoolExecutor(3) as ex:
        res = {n: (st, out) for n, st, out in ex.map(job, JOBS.items())}
    print(f"{'run':32s} {'max|dH/H|':>10s} {'TT':>9s} {'EE':>9s} {'P(k)':>9s}   (vs EFTflag=0, same cosmology)")
    for n, (st, out) in res.items():
        if n.startswith("GR"):
            continue
        if st["status"] != "ok":
            print(f"{n:32s} {st['status']}")
            continue
        # same settings and cosmology (the K-mouflage inis set omnuh2 = 0; ADE runs with the test-ini CMB settings)
        ref = res["GR_nu0" if n.startswith("K") else "GR" if n.startswith("ADE") else "GR_fast"][1]
        r = harness.compare(out, ref, lmax=1900)
        u = r["cl"]["unlensed_scalar"]
        print(f"{n:32s} {r['bg_H']:10.1e} {u['TT']:9.1e} {u['EE']:9.1e} {r['worst_pk']:9.1e}")
    print("expected: all ~1e-6 (coupling 1e-8 / lambda^2 ~ 1e-5); observed: ~1e-4 with the default number of points")


if __name__ == "__main__":
    main()

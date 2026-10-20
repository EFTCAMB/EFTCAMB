"""Verifier (ic1-pidot-transient-noisy-coefficients): independent look at the EFT_IC_type = 1 ingredients.

Uses work/verify/ic1noisy/A0_icdump = A0 + print-only instrumentation in 09_EFTCAMB_IC.f90 case(1) (active only with
EFT_V_ICDUMP=1). For every k at the turn-on it prints C, D1, D2, E, E', the Ridders derivatives C' and D' exactly as the
code computes them (dfridr, h = 0.03 a; D' includes the D1/D2 typo), the Ridders derivative of D1 alone, wide (+-10 %)
centred differences of C and D1, and the pi, pi' actually set. With EFT_V_ICPROF=1 it also prints C(a), D1(a) on
a (1 + j 1e-4), |j| <= 400, through the same EFTpiCfunction / EFTpiDfunction the IC differentiates.

From these we compute d ln pi / d ln a at the turn-on implied by the code's pi' (n_code = pi'/(aH pi), aH = adotoa) and
by a smooth (+-10 %) derivative of the same attractor (n_smooth), and look for discontinuities of C(a) at the background
interpolation nodes.
Usage: python ic1noisy_icdump.py [case ...]
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy"
TREE = W / "A0_icdump"
OUT = W / "icdump"
ENV = {"EFT_V_ICDUMP": "1", "EFT_V_ICPROF": "1"}
CASES = {
    "km_on1e-8": ("5_Kmimic_1", {"EFTCAMB_turn_on_time": "1e-8"}),
    "km_on1.01e-8": ("5_Kmimic_1", {"EFTCAMB_turn_on_time": "1.01e-8"}),
    "km_on1e-9": ("5_Kmimic_1", {"EFTCAMB_turn_on_time": "1e-9"}),
    "km_on1e-7": ("5_Kmimic_1", {"EFTCAMB_turn_on_time": "1e-7"}),
    "km_on1e-9_ai1e-11": ("5_Kmimic_1", {"EFTCAMB_turn_on_time": "1e-9", "model_background_a_ini": "1e-11"}),
    "tg2": ("5_hdsk_TG_2", {}),
    "tg2_on1.01e-8": ("5_hdsk_TG_2", {"EFTCAMB_turn_on_time": "1.01e-8"}),
}
COLS = ["a", "k", "C", "D1", "D2", "E", "Edot", "cd", "dd", "dd1", "cw", "d1w", "Cfun", "D1fun", "pi", "pidot",
        "adotoa", "h0"]


def run(name):
    model, ov = CASES[name]
    out = OUT / f"{name}.npz"
    log = out.with_suffix(".log")
    if not out.exists():
        st = harness.run_model(TREE, harness.PARAMS_DIR / f"{model}.ini", out, overrides=ov, threads=1, timeout=900,
                               extra_env=ENV)
        if st["status"] != "ok":
            return name, st["status"]
    rows, prof = [], []
    for line in log.read_text(errors="replace").splitlines():
        if line.startswith("ICDUMP "):
            rows.append([float(x) for x in line.split()[1:]])
        elif line.startswith("ICPROF "):
            prof.append([float(x) for x in line.split()[1:]])
    return name, (np.array(rows), np.array(prof))


def analyse(name, rows, prof):
    d = {c: rows[:, i] for i, c in enumerate(COLS)}
    order = np.argsort(d["k"])
    d = {c: v[order] for c, v in d.items()}
    k, a, h0, aH = d["k"], d["a"], d["h0"], d["adotoa"]
    den = d["C"] + k**2 * d["D1"] + k**4 * d["D2"]
    pi_chk = -h0 * d["E"] / den
    pidot_chk = h0 * (d["E"] / den**2 * a * aH * (d["cd"] + k**2 * d["dd"]) - d["Edot"] / den)
    pidot_fixtypo = h0 * (d["E"] / den**2 * a * aH * (d["cd"] + k**2 * d["dd1"]) - d["Edot"] / den)
    pidot_smooth = h0 * (d["E"] / den**2 * a * aH * (d["cw"] + k**2 * d["d1w"]) - d["Edot"] / den)
    n_code = d["pidot"] / (aH * d["pi"])
    n_fix = pidot_fixtypo / (aH * pi_chk)
    n_smooth = pidot_smooth / (aH * pi_chk)
    res = {
        "a_on": float(a[0]),
        "n_k": int(len(k)),
        "max_D2_over_D1": float(np.max(np.abs(d["D2"]) / np.maximum(np.abs(d["D1"]), 1e-300))),
        "recon_pi_maxrel": float(np.max(np.abs(pi_chk / d["pi"] - 1))),
        "recon_pidot_maxrel": float(np.max(np.abs(pidot_chk / d["pidot"] - 1))),
        "C_at_on": float(d["C"][0]),
        "C_IC_vs_Cfun_maxrel": float(np.max(np.abs(d["Cfun"] / d["C"] - 1))),
        "D1": float(d["D1"][0]),
        "dlnC_dlna_ridders": float(d["cd"][0] * a[0] / d["C"][0]),
        "dlnC_dlna_wide": float(d["cw"][0] * a[0] / d["C"][0]),
        "dlnD1_dlna_ridders": float(d["dd1"][0] * a[0] / d["D1"][0]),
        "dlnD1_dlna_wide": float(d["d1w"][0] * a[0] / d["D1"][0]),
        "sign_change_den_over_k": int(np.sum(np.diff(np.sign(den)) != 0)),
        "k_sign_change": [float(x) for x in k[1:][np.diff(np.sign(den)) != 0]],
        "median_abs_n_code": float(np.median(np.abs(n_code))),
        "max_abs_n_code": float(np.max(np.abs(n_code))),
        "median_abs_n_fixtypo": float(np.median(np.abs(n_fix))),
        "max_abs_n_fixtypo": float(np.max(np.abs(n_fix))),
        "median_abs_n_smooth": float(np.median(np.abs(n_smooth))),
        "max_abs_n_smooth": float(np.max(np.abs(n_smooth))),
    }
    sel = [int(np.argmin(np.abs(k - kk))) for kk in (0.01, 0.1, 1.0, 5.0, 10.0)]
    res["per_k"] = [{"k": float(k[i]), "den": float(den[i]), "pi": float(d["pi"][i]), "pidot": float(d["pidot"][i]),
                     "n_code": float(n_code[i]), "n_fixtypo": float(n_fix[i]), "n_smooth": float(n_smooth[i])}
                    for i in sel]
    if len(prof):
        x, c, d1 = np.log(prof[:, 0]), prof[:, 1], prof[:, 2]
        # local log-slope from neighbouring points, and its jumps
        sc = np.diff(np.log(np.abs(c))) / np.diff(x)
        sd = np.diff(np.log(np.abs(d1))) / np.diff(x)
        res["prof"] = {
            "C_sign_changes": int(np.sum(np.diff(np.sign(c)) != 0)),
            "dlnC_dlna_local_min_max": [float(sc.min()), float(sc.max())],
            "dlnC_dlna_local_median": float(np.median(sc)),
            "dlnD1_dlna_local_min_max": [float(sd.min()), float(sd.max())],
            "dlnD1_dlna_local_median": float(np.median(sd)),
            "n_local_slope_outliers_C": int(np.sum(np.abs(sc - np.median(sc)) > 1.0)),
            "a_of_C_outliers": [float(v) for v in np.exp(x[1:][np.abs(sc - np.median(sc)) > 1.0])][:12],
        }
    return res


def main():
    names = sys.argv[1:] or list(CASES)
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(run, names))
    summary = {}
    for name in names:
        r = res[name]
        if isinstance(r, str):
            print(name, r)
            continue
        summary[name] = analyse(name, *r)
        print(f"=== {name}")
        print(json.dumps(summary[name], indent=1))
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

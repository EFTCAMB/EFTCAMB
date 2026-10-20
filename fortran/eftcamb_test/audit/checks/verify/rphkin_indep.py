"""Independent check (verifier) of finding rph-kinetic-early-roundoff-ghost.

Claim: for 3_RPH_MB_OmegaDE_1 (alpha_M = Omega_DE, alpha_B = -Omega_DE, alpha_K = alpha_T = 0) the scalar no-ghost
quantity EFT_kinetic (06_abstract_EFTCAMB_model.f90:690) is a difference of two huge terms (4 c M and 8 gamma1 M a^2 h0^2,
RPH builds gamma1 from c) and its sign at a < ~7e-6 is round-off, so the default EFTCAMB_stability_time = 1e-10 rejects a
healthy model.

This script uses code paths different from the repro:
  part "verdict": the real EFTCAMB_Stability_Check via camb.get_background, feedback_level 3, parsing the printed reason
                  (which condition fails, the kinetic value printed, the a where it fails), scanning stability_time, with
                  EFT_ghost_stability on/off, on A0, the t6 fixkin tree (bracket fixed) and B1 (pre-rebase EFTCAMB).
  part "fields":  EFT functions from get_eft_functions (EFTCAMB_eftfunctions, not EFTTestStability); EFT_kinetic is
                  re-evaluated (i) in float64 with the Fortran expression, (ii) in 50-digit mpmath with the code's rounded
                  gamma1, (iii) in mpmath with gamma1 rebuilt exactly from c as RPH defines it, (iv) alpha-basis form
                  18 M^3 H^2 (alpha_K + 6 alpha_B^2).

Run: /Users/marco/Software/miniforge3/bin/python rphkin_indep.py [verdict|fields|all]
"""

import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ROOT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45"
            "/scratchpad/audit")
OUT = ROOT / "work/verify/rphkin"
OUT.mkdir(parents=True, exist_ok=True)
PAR = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")
CHECKS = Path(__file__).resolve().parents[1]
TREES = {"A0": ROOT / "trees/A0", "fixkin": ROOT / "work/t6-stability/fixkin", "B1": ROOT / "trees/B1"}

WORKER = r"""
import sys, json
tree, ini, mode, checks = sys.argv[1:5]
sys.path.insert(0, tree)
sys.path.insert(1, checks)
import numpy as np
import camb
import harness

def make_cp():
    cp = camb.read_ini(ini, no_validate=True)
    keys = harness.parse_ini(ini)
    if cp.EFTCAMB.EFTflag == 0 and int(keys.get("EFTflag", "0")) != 0:
        cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)
    return cp

if mode == "verdict":
    settings = json.loads(sys.argv[5])
    for st, ghost in settings:
        cp = make_cp()
        cp.EFTCAMB.EFTCAMB_stability_time = st
        cp.EFTCAMB.EFT_ghost_stability = ghost
        cp.EFTCAMB.EFTCAMB_feedback_level = 3
        print(f"@@BEGIN {st} {ghost}", flush=True)
        try:
            camb.get_background(cp)
            v = "stable"
        except Exception as e:
            v = "unstable:" + repr(e)[:80]
        print(f"@@END {v}", flush=True)
else:
    a = np.logspace(-10, 0, 1001)
    cp = make_cp()
    cp.EFTCAMB.EFTCAMB_skip_stability = True
    res = camb.get_background(cp)
    fields, arr = cp.EFTCAMB.get_eft_functions(res, a)
    keep = ["adotoa", "Hdot", "EFTc", "EFTOmegaV", "EFTOmegaP", "EFTGamma1V", "EFTGamma2V", "EFTGamma4V", "alphaK",
            "alphaB", "alphaM", "Meff2", "grhov_t", "grhoa2", "EFT_kinetic"]
    out = {k: np.ascontiguousarray(arr[k]) for k in keep}
    out["a"] = a
    out["h0_Mpc"] = np.array(res.Params.EFTCAMB_parameter_cache.h0_Mpc)
    np.savez(sys.argv[5], **out)
    print("@@FIELDS ok", flush=True)
"""


def run_worker(tree, ini, mode, extra):
    env = dict(os.environ, OMP_NUM_THREADS="1", GFORTRAN_UNBUFFERED_PRECONNECTED="y")
    env.pop("PYTHONPATH", None)
    with tempfile.TemporaryDirectory() as d, harness.slot(1):
        p = subprocess.run([sys.executable, "-c", WORKER, str(tree), str(ini), mode, str(CHECKS), extra],
                           capture_output=True, text=True, cwd=d, env=env, timeout=1500,
                           preexec_fn=lambda: os.nice(10))
    return p


def make_ini(model, d):
    ini = Path(d) / f"{model}.ini"
    ini.write_text(f"DEFAULT({PAR}/{model}.ini)\nfeedback_level = 0\nEFTCAMB_write_background = F\n"
                   f"output_root = {d}/out\n")
    return ini


def part_verdict(model="3_RPH_MB_OmegaDE_1", trees=("A0", "fixkin", "B1")):
    sts = [0.1, 1e-3, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 1e-8, 1e-10]
    settings = [[s, True] for s in sts] + [[1e-10, False]]
    result = {}
    with tempfile.TemporaryDirectory() as d:
        ini = make_ini(model, d)
        for t in trees:
            p = run_worker(TREES[t], ini, "verdict", json.dumps(settings))
            txt = p.stdout
            (OUT / f"verdict_{model}_{t}.log").write_text(txt + "\n---STDERR---\n" + p.stderr[-5000:])
            blocks = re.findall(r"@@BEGIN (\S+) (\S+)\n(.*?)@@END (\S+)", txt, re.S)
            rows = []
            for st, ghost, body, v in blocks:
                reasons = sorted(set(re.findall(r"(?:Physical|Mathematical) instability: ([a-z ]+)", body)))
                kin = re.findall(r"ghost instability\. Kinetic term =\s*(\S+)", body)
                where = re.findall(r"Instability detected at a =\s*(\S+)", body)
                rows.append({"stability_time": float(st), "ghost": ghost, "verdict": v.split(":")[0],
                             "reasons": reasons, "kinetic_printed": kin[:3], "a_detected": where[:3]})
            if not blocks:
                rows = {"error": p.stderr[-2000:]}
            result[t] = rows
    return result


def fortran_kin(a, h0, c, om, omp, g1, g2, g4, H, fixed):
    """float64 evaluation of the Fortran expression (A0 form, or the t6 fixkin bracket)."""
    m = 1.0 + om - g4
    if fixed:
        inner = a**2 * h0 * (h0 * (3 * g2**2 + 8 * g1 * m) + 6 * H * g2 * omp)
    else:
        inner = a**2 * h0 * (h0 * (3 * g2**2 + 8 * g1 * m + 6 * H * g2 * omp))
    return 9.0 * m * (4 * c * m + 3 * H**2 * omp**2 * a**2 + inner)


def part_fields(model="3_RPH_MB_OmegaDE_1", tree="A0"):
    import mpmath as mp
    mp.mp.dps = 50
    npz = OUT / f"fields_{model}_{tree}.npz"
    with tempfile.TemporaryDirectory() as d:
        ini = make_ini(model, d)
        p = run_worker(TREES[tree], ini, "fields", str(npz))
        if "@@FIELDS ok" not in p.stdout:
            return {"error": p.stderr[-3000:]}
    f = np.load(npz)
    a, h0 = f["a"], float(f["h0_Mpc"])
    c, om, omp, g1, g2, g4, H = (f[k] for k in ("EFTc", "EFTOmegaV", "EFTOmegaP", "EFTGamma1V", "EFTGamma2V",
                                                "EFTGamma4V", "adotoa"))
    aK, aB, M = f["alphaK"], f["alphaB"], f["Meff2"]
    k64_A0 = fortran_kin(a, h0, c, om, omp, g1, g2, g4, H, fixed=False)
    k64_fix = fortran_kin(a, h0, c, om, omp, g1, g2, g4, H, fixed=True)
    ana = 18 * M**3 * H**2 * (aK + 6 * aB**2)
    big = 36 * M**2 * np.abs(c)
    rows = []
    kmp_code_g1, kmp_exact_A0, kmp_exact_fix = [], [], []
    for i in range(len(a)):
        A, h, C, O, OP, G1, G2, G4, HH, AK = (mp.mpf(float(x)) for x in
                                              (a[i], h0, c[i], om[i], omp[i], g1[i], g2[i], g4[i], H[i], aK[i]))
        m = 1 + O - G4
        def kin(gam1, fixed):
            if fixed:
                inner = A**2 * h * (h * (3 * G2**2 + 8 * gam1 * m) + 6 * HH * G2 * OP)
            else:
                inner = A**2 * h * (h * (3 * G2**2 + 8 * gam1 * m + 6 * HH * G2 * OP))
            return 9 * m * (4 * C * m + 3 * HH**2 * OP**2 * A**2 + inner)
        # RPH (007p2_RPH.f90:1749): gamma1 = (alpha_K (1+PM) H^2 - 2c) / (4 h0^2 a^2); here alpha_T = 0 -> PM = OmegaV
        g1_exact = (AK * (1 + O) * HH**2 - 2 * C) / (4 * h**2 * A**2)
        kmp_code_g1.append(float(kin(G1, False)))
        kmp_exact_A0.append(float(kin(g1_exact, False)))
        kmp_exact_fix.append(float(kin(g1_exact, True)))
    kmp_code_g1, kmp_exact_A0, kmp_exact_fix = map(np.array, (kmp_code_g1, kmp_exact_A0, kmp_exact_fix))
    r = ana / big
    neg64 = k64_A0 < 0
    negfix = k64_fix < 0
    out = {
        "tree": tree, "n": int(len(a)), "h0_Mpc": h0,
        "get_eft_functions_EFT_kinetic_all_zero": bool(np.all(f["EFT_kinetic"] == 0)),
        "n_neg_float64_A0form": int(neg64.sum()), "n_neg_float64_fixform": int(negfix.sum()),
        "a_range_neg_float64_A0form": [float(a[neg64].min()), float(a[neg64].max())] if neg64.any() else None,
        "a_range_neg_float64_fixform": [float(a[negfix].min()), float(a[negfix].max())] if negfix.any() else None,
        "n_neg_mp_code_gamma1": int((kmp_code_g1 < 0).sum()),
        "n_neg_mp_exact_gamma1_A0form": int((kmp_exact_A0 < 0).sum()),
        "n_neg_mp_exact_gamma1_fixform": int((kmp_exact_fix < 0).sum()),
        "n_neg_analytic": int((ana < 0).sum()),
        "max_ratio_true_over_cancelling_at_neg_A0form": float(r[neg64].max()) if neg64.any() else None,
        "min_ratio_true_over_cancelling_at_nonneg": float(r[~neg64].min()),
        "max_abs_k64_over_cancelling_at_neg": float((np.abs(k64_A0) / big)[neg64].max()) if neg64.any() else None,
        "exact_fixform_over_analytic_minmax": [float(np.min(kmp_exact_fix / ana)), float(np.max(kmp_exact_fix / ana))],
        "exact_A0form_over_analytic_minmax": [float(np.min(kmp_exact_A0 / ana)), float(np.max(kmp_exact_A0 / ana))],
        "a_where_true_over_cancelling_eq_1e-15": float(np.interp(-15, np.log10(r[(a > 1e-8)]), a[(a > 1e-8)])),
    }
    for i in (0, 100, 150, 200, 250, 300, 350, 400, 450, 470, 500, 530, 600, 700, 800, 1000):
        rows.append({"a": float(a[i]), "k64_A0form": float(k64_A0[i]), "k64_fixform": float(k64_fix[i]),
                     "mp_code_gamma1": float(kmp_code_g1[i]), "mp_exact_gamma1_A0form": float(kmp_exact_A0[i]),
                     "mp_exact_gamma1_fixform": float(kmp_exact_fix[i]), "analytic": float(ana[i]),
                     "cancelling_36M2c": float(big[i]), "true_over_cancelling": float(r[i]),
                     "OmegaDE_proxy_alphaB": float(aB[i]), "alphaM": float(f["alphaM"][i])})
    out["rows"] = rows
    return out


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    res = {}
    if what in ("fields", "all"):
        res["fields_A0"] = part_fields(tree="A0")
    if what in ("verdict", "all"):
        res["verdict"] = part_verdict()
    (OUT / f"indep_{what}.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))

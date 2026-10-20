"""Verify t3-bg-table-resolution-gr-floor at the C_l / P(k) level with independent comparisons.

1. self-convergence in model_background_num_points (no GR reference): K-mouflage eps2_0=-1e-8, N=1000,2000,4000 vs
   N=10000 -> expect ~N^-2 scaling.
2. coupling independence: eps2_0 = -1e-6, -1e-8, -1e-10 at default N vs EFTflag=0 with MATCHED accuracy settings
   (the repro uses IntTolBoost=2 for the GR reference only).
3. AccuracyBoost=2 on both model and GR: does the floor go away? (num_points is not tied to AccuracyBoost.)
4. CAMB's own default numerical error for context: GR default vs GR AccuracyBoost=2,lAccuracyBoost=2, and
   K-mouflage default (N=1000) vs K-mouflage high accuracy (N=10000, AccuracyBoost=2, lAccuracyBoost=2).
5. quintessence lambda scan at default N (floor) and N=11000.
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

P = harness.PARAMS_DIR
OUT = harness.AUDIT_ROOT / "work" / "verify" / "bgtable_cl"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "2000", "num_redshiftwindows": "0"}
K = P / "5_Kmouflage_1.ini"
Q = P / "5_quint_exponential.ini"
G = P / "1_EFT_GR.ini"
HI = {"Accuracy.AccuracyBoost": 2.0, "Accuracy.lAccuracyBoost": 2.0}
AB2 = {"Accuracy.AccuracyBoost": 2.0}

JOBS = {
    "k_e8_N1000": (K, {"eps2_0": "-1e-08"}, {}),
    "k_e8_N2000": (K, {"eps2_0": "-1e-08", "model_background_num_points": "2000"}, {}),
    "k_e8_N4000": (K, {"eps2_0": "-1e-08", "model_background_num_points": "4000"}, {}),
    "k_e8_N10000": (K, {"eps2_0": "-1e-08", "model_background_num_points": "10000"}, {}),
    "k_e6_N1000": (K, {"eps2_0": "-1e-06"}, {}),
    "k_e10_N1000": (K, {"eps2_0": "-1e-10"}, {}),
    "gr_nu0": (G, {"omnuh2": "0.0"}, {}),
    "gr_nu0_itb2": (G, {"omnuh2": "0.0"}, {"Accuracy.IntTolBoost": 2.0}),
    "k_e8_N1000_ab2": (K, {"eps2_0": "-1e-08"}, AB2),
    "gr_nu0_ab2": (G, {"omnuh2": "0.0"}, AB2),
    "k_e8_N10000_hi": (K, {"eps2_0": "-1e-08", "model_background_num_points": "10000"}, HI),
    "gr_nu0_hi": (G, {"omnuh2": "0.0"}, HI),
    "q_l3e-3_N1100": (Q, {"EFT_lambda": "0.003"}, {}),
    "q_l1e-2_N1100": (Q, {"EFT_lambda": "0.01"}, {}),
    "q_l1e-3_N1100": (Q, {"EFT_lambda": "0.001"}, {}),
    "q_l3e-3_N11000": (Q, {"EFT_lambda": "0.003", "model_background_num_points": "11000"}, {}),
    "q_l1e-3_N11000": (Q, {"EFT_lambda": "0.001", "model_background_num_points": "11000"}, {}),
    "gr": (G, {}, {}),
}

PAIRS = [
    ("self-conv", "k_e8_N1000", "k_e8_N10000"),
    ("self-conv", "k_e8_N2000", "k_e8_N10000"),
    ("self-conv", "k_e8_N4000", "k_e8_N10000"),
    ("coupling", "k_e6_N1000", "gr_nu0"),
    ("coupling", "k_e8_N1000", "gr_nu0"),
    ("coupling", "k_e10_N1000", "gr_nu0"),
    ("coupling", "k_e8_N10000", "gr_nu0"),
    ("ref-ITB", "k_e8_N10000", "gr_nu0_itb2"),
    ("ref-ITB", "gr_nu0", "gr_nu0_itb2"),
    ("AB2", "k_e8_N1000_ab2", "gr_nu0_ab2"),
    ("camb-err", "gr_nu0", "gr_nu0_hi"),
    ("camb-err", "gr_nu0_ab2", "gr_nu0_hi"),
    ("total-err", "k_e8_N1000", "k_e8_N10000_hi"),
    ("total-err", "k_e8_N10000_hi", "gr_nu0_hi"),
    ("quint", "q_l1e-2_N1100", "gr"),
    ("quint", "q_l3e-3_N1100", "gr"),
    ("quint", "q_l1e-3_N1100", "gr"),
    ("quint", "q_l3e-3_N11000", "gr"),
    ("quint", "q_l1e-3_N11000", "gr"),
    ("quint-self", "q_l3e-3_N1100", "q_l3e-3_N11000"),
]


def job(item):
    name, (ini, ov, ps) = item
    out = OUT / f"{name}.npz"
    if out.exists() and json.loads(out.with_suffix(".json").read_text())["status"] == "ok":
        return name, {"status": "ok"}
    st = harness.run_model("A0", ini, out, overrides={**FAST, **ov}, python_set=ps, threads=1, timeout=1500)
    return name, st


def main() -> None:
    sel = sys.argv[1:] or list(JOBS)
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(job, [(n, JOBS[n]) for n in sel]))
    for n, st in res.items():
        if st["status"] != "ok":
            print(n, st["status"], st.get("log_tail", [])[-5:])
    rows = []
    print(f"{'group':10s} {'A':18s} {'B':16s} {'dH/H':>9s} {'TT':>9s} {'EE':>9s} {'TE':>9s} {'P(k)':>9s} {'s8':>9s}")
    for g, a, b in PAIRS:
        fa, fb = OUT / f"{a}.npz", OUT / f"{b}.npz"
        if not (fa.exists() and fb.exists()):
            continue
        r = harness.compare(fa, fb, lmax=1900)
        u = r["cl"]["unlensed_scalar"]
        row = dict(group=g, a=a, b=b, bg_H=r["bg_H"], TT=u["TT"], EE=u["EE"], TE=u["TE"], pk=r["worst_pk"],
                   sigma8=r.get("sigma8", float("nan")))
        rows.append(row)
        print(f"{g:10s} {a:18s} {b:16s} {row['bg_H']:9.2e} {row['TT']:9.2e} {row['EE']:9.2e} {row['TE']:9.2e} "
              f"{row['pk']:9.2e} {row['sigma8']:9.2e}")
    (OUT / "summary.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()

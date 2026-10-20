"""Verifier (spline-integral-convention): model-level check through the full CAMB/EFTCAMB pipeline (harness).

Consumers of parametrized_function_1D%integral are only designer f(R) (007p3) and shift-symmetric alpha_B (007p7).
Constant w = -0.9 written as a spline w(a) must reproduce the EFTwDE = 1 (constant w0 = -0.9) run.
Cases (tree A0, release):
  DFR_ref     4_Des_F_R_wCDM_1 as is (EFTwDE=1, w0=-0.9)
  DFR_s9      EFTwDE=9 nodes 1e-8, 0.2..1        (reporter's 2nd layout)
  DFR_s9b     EFTwDE=9 nodes 1e-9, 1e-6, 1e-3, 0.3, 0.7, 1.2 (last node beyond a = 1, first below a_ini)
  DFR_s9def   EFTwDE=9 default nodes (x_i = (i-1)/(N-1), x1 = 0)
  DFR_s10     EFTwDE=10 Logx nodes 1e-8..1
  DFR_s10b    EFTwDE=10 Logx nodes 1e-9 .. 1.2 (last node beyond a = 1)
  SS_ref      007p7 shift-symmetric alpha_B with EFTwDE=1, w0=-0.9
  SS_s9b      same with EFTwDE=9 nodes as DFR_s9b
  SS_s10b     same with EFTwDE=10 nodes as DFR_s10b
For runs that complete, compares eft_grhov_t (= grhov*integral(a) in both models) with the reference and spectra.
usage: /Users/marco/Software/miniforge3/bin/python verify_spline_integral_models.py
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "spline_integral" / "models"
DFR = harness.PARAMS_DIR / "4_Des_F_R_wCDM_1.ini"
SS = harness.REPO / "fortran" / "eftcamb_test" / "audit" / "inis" / "007p7-shiftsym.ini"
W = -0.9


def spline(flag, xs, logx=False, nodes=True):
    ov = {"EFTwDE": str(flag), "EFTw_Spline_Pixels": str(len(xs)), "EFTw_null_value": repr(W)}
    if logx:
        ov["EFTw_Logx_Flag"] = "T"
    for i, x in enumerate(xs):
        if nodes:
            ov[f"EFTwx{i + 1}"] = repr(float(x))
        ov[f"EFTwv{i + 1}"] = repr(W)
    return ov


N9 = [1e-8, 0.2, 0.4, 0.6, 0.8, 1.0]
N9B = [1e-9, 1e-6, 1e-3, 0.3, 0.7, 1.2]
N10 = [10.0 ** e for e in range(-8, 1)]
N10B = [1e-9, 1e-7, 1e-5, 1e-3, 1e-2, 0.1, 0.4, 0.8, 1.2]
CASES = {
    "DFR_ref": (DFR, {}),
    "DFR_s9": (DFR, spline(9, N9)),
    "DFR_s9b": (DFR, spline(9, N9B)),
    "DFR_s9def": (DFR, spline(9, [0, 0.2, 0.4, 0.6, 0.8, 1.0], nodes=False)),
    "DFR_s10": (DFR, spline(10, N10, logx=True)),
    "DFR_s10b": (DFR, spline(10, N10B, logx=True)),
    "SS_ref": (SS, {"EFTwDE": "1", "EFTw0": repr(W)}),
    "SS_s9b": (SS, spline(9, N9B)),
    "SS_s10b": (SS, spline(10, N10B, logx=True)),
}


def run(name):
    ini, ov = CASES[name]
    out = OUT / f"{name}.npz"
    js = out.with_suffix(".json")
    if js.exists():
        return name, json.loads(js.read_text())
    return name, harness.run_model("A0", ini, out, overrides=ov, threads=1, timeout=900)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(4) as ex:
        sts = dict(ex.map(run, CASES))
    for name, st in sts.items():
        tail = ""
        if st["status"] != "ok":
            errs = [ln for ln in st.get("log_tail", []) if "rror" in ln or "FAIL" in ln or "NaN" in ln]
            tail = " | " + " / ".join(e.strip()[:150] for e in errs[-3:])
        print(f"{name:10s} {st['status']}{tail}")
    print()
    for ref, tests in (("DFR_ref", ["DFR_s9", "DFR_s9b", "DFR_s9def", "DFR_s10", "DFR_s10b"]),
                       ("SS_ref", ["SS_s9b", "SS_s10b"])):
        if sts[ref]["status"] != "ok":
            print(ref, "reference failed")
            continue
        r = harness.load(OUT / f"{ref}.npz")
        a = r["eft_a_grid"]
        for t in tests:
            if sts[t]["status"] != "ok":
                continue
            d = harness.load(OUT / f"{t}.npz")
            c = harness.compare(d, r)
            g, g0 = d["eft_grhov_t"], r["eft_grhov_t"]
            m = (a > 1e-6) & (np.abs(g0) > 0)
            rel = np.abs(g[m] - g0[m]) / np.abs(g0[m])
            k = np.argmax(rel)
            print(f"{t:10s} vs {ref}: max rel |grhov_t| diff {rel[k]:.3e} at a={a[m][k]:.3e} "
                  f"(test {g[m][k]:.4e} ref {g0[m][k]:.4e}); grhov_t(a=1) test {g[-1]:.5e} ref {g0[-1]:.5e}; "
                  f"worst Cl {c['worst_cl']:.3e} P(k) {c['worst_pk']:.3e} sigma8 {c['sigma8']:.3e}")


if __name__ == "__main__":
    main()

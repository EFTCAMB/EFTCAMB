"""End-to-end check of finding pade-fourth-derivative through the full Pure-EFT model (harness, A0 vs one-line fix).

Pure EFT (EFTflag=1, PureEFTmodel=1) on LCDM background with Omega(a) = Pade [N=2/D=1] around a0=0:
    Omega(a) = (0.05 a + 0.03 a^2) / (1 + 0.5 a)
Compares eft_EFTOmegaPPPP (via get_eft_functions) to the exact sympy derivative and the A0 / fixed runs.
"""

import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify")
FIX = W / "pade_fix"
INI = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"
OUT = W / "pade_e2e"

A0, UP, DOWN = 0.0, [0.0, 0.05, 0.03], [0.5]
OV = {"PureEFTmodelOmega": "6", "EFTOmega_Pade_order_N": str(len(UP) - 1), "EFTOmega_Pade_order_D": str(len(DOWN)),
      "EFTOmegaa0": str(A0)}
OV.update({f"EFTOmegaN{i}": str(v) for i, v in enumerate(UP)})
OV.update({f"EFTOmegaD{i + 1}": str(v) for i, v in enumerate(DOWN)})


def exact():
    a = sp.symbols("a")
    F = sum(sp.Rational(str(c)) * (a - A0) ** i for i, c in enumerate(UP)) / (
        1 + sum(sp.Rational(str(c)) * (a - A0) ** (i + 1) for i, c in enumerate(DOWN)))
    return [sp.lambdify(a, sp.diff(F, a, k), "numpy") for k in range(5)]


def main():
    from concurrent.futures import ThreadPoolExecutor
    jobs = {"A0": "A0", "FIX": str(FIX)}
    with ThreadPoolExecutor(2) as ex:
        st = dict(zip(jobs, ex.map(lambda t: harness.run_model(t, INI, OUT / f"{[k for k, v in jobs.items() if v == t][0]}.npz",
                                                               overrides=OV, threads=1), jobs.values())))
    for k, s in st.items():
        print(k, s["status"], s.get("log_tail", "")[-3:] if s["status"] != "ok" else "")
    d = {k: harness.load(OUT / f"{k}.npz") for k in jobs}
    ex_f = exact()
    ag = d["A0"]["eft_a_grid"]
    m = ag > 1e-3
    res = {}
    for name, key, k in [("Omega", "eft_EFTOmegaV", 0), ("OmegaP", "eft_EFTOmegaP", 1), ("OmegaPP", "eft_EFTOmegaPP", 2),
                         ("OmegaPPP", "eft_EFTOmegaPPP", 3), ("OmegaPPPP", "eft_EFTOmegaPPPP", 4)]:
        e = ex_f[k](ag) * np.ones_like(ag)
        for t in jobs:
            v = d[t][key]
            r = float(np.max(np.abs(v[m] - e[m])) / np.max(np.abs(e[m])))
            res[f"{t}_{name}_maxrel_vs_exact"] = r
    for t in jobs:
        i1 = np.argmin(abs(ag - 1.0))
        res[f"{t}_OmegaPPPP_at_a1"] = float(d[t]["eft_EFTOmegaPPPP"][i1])
    res["exact_OmegaPPPP_at_a1"] = float(ex_f[4](1.0))
    for f in ["EFTcdotdot", "EFTLambdadotdot", "EFTc", "EFTLambda", "EFTcdot", "EFTLambdadot"]:
        a_, b_ = d["A0"][f"eft_{f}"], d["FIX"][f"eft_{f}"]
        res[f"A0_vs_FIX_{f}"] = float(np.max(np.abs(a_ - b_)) / max(np.max(np.abs(b_)), 1e-300))
    c = harness.compare(d["A0"], d["FIX"])
    res["A0_vs_FIX_worst_cl"] = c["worst_cl"]
    res["A0_vs_FIX_worst_pk"] = c["worst_pk"]
    res["A0_vs_FIX_bit_identical"] = c["bit_identical"]
    print(json.dumps(res, indent=1, default=str))
    (OUT / "summary.json").write_text(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()

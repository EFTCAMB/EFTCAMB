"""Track t4a-param: before/after model-level numbers for the parametrization findings (A0 vs the FIX tree copy built
with make_fix_tree.py).  Each case runs on both trees; 'ref' cases are also compared with an exact-equivalent
reference run on A0 (same function written with an analytic parametrization).

usage: python fix_compare.py [--only substr] [--jobs 3] [--report]
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import harness  # noqa: E402
import model_checks as mc  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "t4a-param" / "fixcmp"
FIX = harness.AUDIT_ROOT / "work" / "t4a-param" / "FIX"
PEFT_LIN = mc.PEFT_LIN
NODES_PART = [0.01, 0.1, 0.25, 0.4, 0.6, 0.9]

CASES = {k: mc.CASES[k][:2] for k in ["peft_spline_lin", "peft_spline5_o3_lin", "peft_steplog2", "dfr_steplog_step",
                                      "dfr_spline_const_x1", "dfr_spline5_const"]}
# w = -1 written as splines whose nodes do not cover [1e-8, 1] (null_value = -1): must equal EFTwDE = 0
CASES["peft_w_ref"] = (PEFT_LIN, {"EFTwDE": 0})
CASES["peft_w_spline_part"] = (PEFT_LIN, {"EFTwDE": 9, **mc.spl("EFTw", NODES_PART, [-1.0] * 6,
                                                                {"EFTw_null_value": -1.0})})
CASES["peft_w_spline5_part"] = (PEFT_LIN, {"EFTwDE": 10, **mc.spl("EFTw", NODES_PART, [-1.0] * 6,
                                                                  {"EFTw_null_value": -1.0, "EFTw_Spline_Order": 5})})
# Omega as a genuinely rational Pade function (fourth derivative matters): A0 vs FIX
PADE = {"PureEFTmodelOmega": 6, "EFTOmega_Pade_order_N": 2, "EFTOmega_Pade_order_D": 1, "EFTOmegaa0": 1.0,
        "EFTOmegaN0": 0.05, "EFTOmegaN1": 0.05, "EFTOmegaN2": 0.02, "EFTOmegaD1": 0.3}
CASES["peft_pade_curved"] = (PEFT_LIN, {**PADE, "EFT_mass_stability": "F"})
CASES["peft_pade_curved_mstab"] = (PEFT_LIN, PADE)
REF = {"peft_spline_lin": "model", "peft_spline5_o3_lin": "model", "dfr_spline_const_x1": "model",
       "dfr_spline5_const": "model", "peft_w_spline_part": "peft_w_ref", "peft_w_spline5_part": "peft_w_ref"}


def run(item):
    name, tag = item
    ini, ov = CASES[name]
    tree = "A0" if tag == "A0" else FIX
    out = W / tag / f"{name}.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    st = harness.run_model(tree, ini, out, overrides={k: str(v) for k, v in ov.items()}, threads=2, timeout=1500)
    return name, tag, st["status"]


def short(cmp):
    eft = sorted(cmp["eft"].items(), key=lambda kv: -kv[1])[:3]
    return (f"worst_cl={cmp['worst_cl']:.3e} worst_pk={cmp['worst_pk']:.3e} s8={cmp['sigma8']:.2e} "
            f"eft={[(k, f'{v:.1e}') for k, v in eft]}")


def report():
    rows = {}
    for name in CASES:
        row = {}
        for tag in ("A0", "FIX"):
            j = W / tag / f"{name}.json"
            row[tag] = json.loads(j.read_text())["status"] if j.exists() else "missing"
        line = f"{name:22s} A0={row['A0']:12s} FIX={row['FIX']:12s}"
        p = {t: W / t / f"{name}.npz" for t in ("A0", "FIX")}
        if name in REF:
            ref = mc.ref_npz(CASES[name][0]) if REF[name] == "model" else W / "A0" / f"{REF[name]}.npz"
            for tag in ("A0", "FIX"):
                if row[tag] == "ok" and ref.exists():
                    c = harness.compare(harness.load(p[tag]), harness.load(ref))
                    row[f"{tag}_vs_ref"] = c
                    line += f"\n    {tag} vs exact-equivalent reference: {short(c)}"
        elif row["A0"] == "ok" and row["FIX"] == "ok":
            c = harness.compare(harness.load(p["A0"]), harness.load(p["FIX"]))
            row["A0_vs_FIX"] = c
            line += f"\n    A0 vs FIX: {short(c)}"
        rows[name] = row
        print(line)
    (W / "summary.json").write_text(json.dumps(rows, indent=1, default=str))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()
    if args.report:
        report()
        return
    items = [(n, t) for n in CASES if args.only in n for t in ("A0", "FIX")
             if not (n == "peft_w_ref" and t == "FIX")]
    items = [it for it in items if not (W / it[1] / f"{it[0]}.json").exists()]
    with ThreadPoolExecutor(args.jobs) as ex:
        for name, tag, status in ex.map(run, items):
            print(name, tag, status, flush=True)


if __name__ == "__main__":
    main()

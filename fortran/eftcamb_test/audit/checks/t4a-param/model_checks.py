"""Track t4a-param: model-level consequences of parametrization defects.

Each case re-expresses a function of a reference test model with a different 1D parametrization that represents
the SAME function exactly (constant w, linear Omega); the spectra must then reproduce the reference run.
Cases marked 'fixcmp' are instead compared between the A0 build and a fixed tree copy.

usage: python model_checks.py --tree A0|PATH [--only substr] [--ref-tree A0] [--tag NAME] [--jobs 4]
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

PAR = harness.PARAMS_DIR
WORK = harness.AUDIT_ROOT / "work" / "t4a-param" / "models"

DES_FR = PAR / "4_Des_F_R_wCDM_1.ini"  # EFTwDE = 1, EFTw0 = -0.9
PEFT_LIN = PAR / "2_PEFT_Omega_lin_1.ini"  # Omega = 1.0 * a

X6 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
LOGX9 = [10.0 ** e for e in range(-8, 1)]


def spl(prefix: str, xs, vs, extra=None):
    d = {f"{prefix}_Spline_Pixels": len(xs)}
    d.update({f"{prefix}x{i + 1}": x for i, x in enumerate(xs)})
    d.update({f"{prefix}v{i + 1}": v for i, v in enumerate(vs)})
    d.update(extra or {})
    return d


CASES = {
    # ---- designer f(R): constant w = -0.9 written with the allocator parametrizations (EFTwDE 6..11)
    "dfr_pade_const": (DES_FR, {"EFTwDE": 6, "EFTw_Pade_order_N": 0, "EFTw_Pade_order_D": 0, "EFTwa0": 0.0,
                                "EFTwN0": -0.9}, "ref"),
    "dfr_fourier_const": (DES_FR, {"EFTwDE": 7, "EFTw_Fourier_order_cos": 0, "EFTw_Fourier_order_sin": 0,
                                   "EFTw_a0": -0.9, "EFTw_phase": 0.0}, "ref"),
    "dfr_steplog_const": (DES_FR, {"EFTwDE": 8, "EFTw_v1": -0.9, "EFTw_v2": -0.9, "EFTw_at": 0.5,
                                   "EFTw_delta": 0.3}, "ref"),
    "dfr_spline_const": (DES_FR, {"EFTwDE": 9, **spl("EFTw", X6, [-0.9] * 6, {"EFTw_null_value": -0.9})}, "ref"),
    "dfr_spline_const_x1": (DES_FR, {"EFTwDE": 9, **spl("EFTw", [1e-8] + X6[1:], [-0.9] * 6,
                                                         {"EFTw_null_value": -0.9})}, "ref"),
    "dfr_spline5_const": (DES_FR, {"EFTwDE": 10, **spl("EFTw", LOGX9, [-0.9] * 9,
                                                       {"EFTw_null_value": -0.9, "EFTw_Logx_Flag": "T",
                                                        "EFTw_Spline_Order": 5})}, "ref"),
    "dfr_exp2_const": (DES_FR, {"EFTwDE": 11, "EFTw0": -0.9, "EFTw1": 0.0}, "ref"),
    # steplog with a real transition: A0 vs fixed integral
    "dfr_steplog_step": (DES_FR, {"EFTwDE": 8, "EFTw_v1": -1.0, "EFTw_v2": -0.8, "EFTw_at": 0.3,
                                  "EFTw_delta": 0.5}, "fixcmp"),
    # ---- pure EFT: Omega = a written as splines of exactly linear data
    "peft_spline_lin": (PEFT_LIN, {"PureEFTmodelOmega": 9, **spl("EFTOmega", X6, X6)}, "ref"),
    "peft_spline5_lin": (PEFT_LIN, {"PureEFTmodelOmega": 10, **spl("EFTOmega", X6, X6,
                                                                   {"EFTOmega_Spline_Order": 5})}, "ref"),
    "peft_spline5_o3_lin": (PEFT_LIN, {"PureEFTmodelOmega": 10, **spl("EFTOmega", X6, X6,
                                                                      {"EFTOmega_Spline_Order": 3})}, "ref"),
    "peft_pade_lin": (PEFT_LIN, {"PureEFTmodelOmega": 6, "EFTOmega_Pade_order_N": 1, "EFTOmega_Pade_order_D": 0,
                                 "EFTOmegaa0": 0.0, "EFTOmegaN0": 0.0, "EFTOmegaN1": 1.0}, "ref"),
    # ---- pure EFT: Omega = 0.1 a^1.5 as power law (reference) and as quintic spline in ln a (poly5 bug)
    "peft_pow15_ref": (PEFT_LIN, {"EFT_mass_stability": "F", "PureEFTmodelOmega": 3, "EFTOmega0": 0.1, "EFTOmegaExp": 1.5}, "fixcmp"),
    "peft_spline5_pow15": (PEFT_LIN, {"EFT_mass_stability": "F", "PureEFTmodelOmega": 10, **spl("EFTOmega", [10.0 ** (-4 + 0.5 * i) for i in range(9)],
                                                                    [0.1 * (10.0 ** (-4 + 0.5 * i)) ** 1.5 for i in range(9)],
                                                                    {"EFTOmega_Spline_Order": 5, "EFTOmega_Logx_Flag": "T"})},
                           "fixcmp"),
    "peft_spline5_pow15_lin": (PEFT_LIN, {"EFT_mass_stability": "F", "PureEFTmodelOmega": 10, **spl("EFTOmega", [0.125 * i for i in range(9)],
                                                                        [0.1 * (0.125 * i) ** 1.5 for i in range(9)],
                                                                        {"EFTOmega_Spline_Order": 5})}, "fixcmp"),
    "peft_spline5_pow15_logx": (PEFT_LIN, {"EFT_mass_stability": "F", "PureEFTmodelOmega": 10, **spl("EFTOmega", [10.0 ** (-2 + 0.25 * i) for i in range(9)],
                                                                    [0.1 * (10.0 ** (-2 + 0.25 * i)) ** 1.5 for i in range(9)],
                                                                    {"EFTOmega_Spline_Order": 5, "EFTOmega_Logx_Flag": "T"})},
                                "fixcmp"),
    "peft_steplog2": (PEFT_LIN, {"EFT_mass_stability": "F", "PureEFTmodelOmega": 8, "EFTOmega_v1": 0.0, "EFTOmega_v2": 1.5,
                                 "EFTOmega_at": 0.5, "EFTOmega_delta": 0.5}, "fixcmp"),
    "peft_steplog3": (PEFT_LIN, {"EFT_mass_stability": "F", "PureEFTmodelOmega": 8, "EFTOmega_v1": 0.0, "EFTOmega_v2": 1.0,
                                 "EFTOmega_at": 1.0, "EFTOmega_delta": 1.0}, "fixcmp"),
    # ---- pure EFT: Omega steplog (fourth derivative NaN at a > aT exp(delta^2))
    "peft_steplog": (PEFT_LIN, {"PureEFTmodelOmega": 8, "EFTOmega_v1": 0.0, "EFTOmega_v2": 0.05,
                                "EFTOmega_at": 0.5, "EFTOmega_delta": 0.3}, "fixcmp"),
}


def ref_npz(ini: Path) -> Path:
    return harness.AUDIT_ROOT / "runs" / "A0" / f"{ini.stem}.npz"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default="A0")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--only", default="")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--compare", default="", help="TAG1,TAG2: compare finished runs of two tags")
    args = ap.parse_args()
    if args.compare:
        t1, t2 = args.compare.split(",")
        for name in CASES:
            if args.only not in name:
                continue
            p1, p2 = WORK / t1 / f"{name}.npz", WORK / t2 / f"{name}.npz"
            s1 = json.loads((WORK / t1 / f"{name}.json").read_text())["status"] if (WORK / t1 / f"{name}.json").exists() else "missing"
            s2 = json.loads((WORK / t2 / f"{name}.json").read_text())["status"] if (WORK / t2 / f"{name}.json").exists() else "missing"
            row = {"status": (s1, s2)}
            if p1.exists() and p2.exists():
                cmp = harness.compare(harness.load(p1), harness.load(p2))
                row.update(worst_cl=cmp["worst_cl"], worst_pk=cmp["worst_pk"], sigma8=cmp.get("sigma8"),
                           eft={k: v for k, v in sorted(cmp["eft"].items(), key=lambda kv: -kv[1])[:5]})
            print(name, json.dumps(row, default=str)[:700])
        return
    tag = args.tag or Path(args.tree).name
    outdir = WORK / tag
    outdir.mkdir(parents=True, exist_ok=True)
    sel = [k for k in CASES if args.only in k]

    def one(name):
        ini, ov, kind = CASES[name]
        st = harness.run_model(args.tree, ini, outdir / f"{name}.npz", overrides={k: str(v) for k, v in ov.items()},
                               threads=2, timeout=1500)
        return name, st

    with ThreadPoolExecutor(args.jobs) as ex:
        results = dict(ex.map(one, sel))
    summary = {}
    for name in sel:
        ini, ov, kind = CASES[name]
        st = results[name]
        row = {"status": st["status"], "log_tail": st.get("log_tail", "")[-600:]}
        if st["status"] == "ok" and kind == "ref":
            cmp = harness.compare(harness.load(outdir / f"{name}.npz"), harness.load(ref_npz(ini)))
            row.update(worst_cl=cmp["worst_cl"], worst_pk=cmp["worst_pk"], sigma8=cmp["sigma8"],
                       bg_H=cmp["bg_H"], bit_identical=cmp["bit_identical"],
                       eft={k: v for k, v in sorted(cmp["eft"].items(), key=lambda kv: -kv[1])[:6]})
        summary[name] = row
        print(name, json.dumps(row, default=str)[:900])
    (outdir / "summary.json").write_text(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()

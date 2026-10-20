"""Track 5 (numerical robustness): run every test model with the EFTCAMB control parameters varied.

Variations (relative to the value in each model ini, or the code default when the ini does not set it):
    AB2      accuracy_boost = 2                       (noise floor)
    bto_x10  EFTCAMB_back_turn_on * 10, bto_d10 / 10
    pto_x10  EFTCAMB_turn_on_time * 10, pto_d10 / 10 (perturbation turn-on)
    grt_x10  EFTCAMB_GR_threshold * 10, grt_d10 / 10
    stt_x10  EFTCAMB_stability_time * 10, stt_d10 / 10
    ic2/ic3  EFT_IC_type = 2 / 3
    ref      default settings rerun (determinism check against runs/A0)

Usage: python run_variations.py VAR[,VAR...] [--models pattern] [--jobs 4] [--tree A0]
Outputs: $WORK/runs/<VAR>/<model>.npz (skips existing outputs).
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t5-robust"
DEFAULTS = {
    "EFTCAMB_back_turn_on": 1e-8,
    "EFTCAMB_turn_on_time": 1e-2,
    "EFTCAMB_GR_threshold": 1e-8,
    "EFTCAMB_stability_time": 1e-10,
}
# representative models of every family (used when running all 73 models is too slow)
REPRESENTATIVE = [
    "2_PEFT_Omega_const_2", "2_PEFT_Omega_exp_1", "2_PEFT_Omega_lin_2", "2_PEFT_Omega_power_2", "2_PEFT_gamma1_const_2",
    "2_PEFT_gamma3_lin_2", "2_PEFT_gamma3_power_2", "2_PEFT_gamma4_const_1", "2_PEFT_gamma5_exp_2",
    "2_PEFT_gamma6_lin_1", "2_PEFT_wDE_wCDM_1", "2_PEFT_wDE_LCDM_1", "3_RPH_KT_power_1", "3_RPH_MK_power_1",
    "3_RPH_MK_OmegaDE_1", "4_Des_F_R_LCDM_2", "4_Des_F_R_wCDM_1", "4_Des_mc5e_CPL_1", "4_Des_mc5e_wCDM_1",
    "5_ADE_1", "5_Horava_eta_1", "5_Kmimic_1", "5_Kmouflage_3", "5_quint_axion", "5_quint_exponential", "5_scg",
    "5_hdsk_TG_1", "5_hdsk_TG_2",
]
# smaller core set used when the shared machine is overloaded
CORE = [
    "2_PEFT_Omega_const_2", "2_PEFT_Omega_exp_1", "2_PEFT_Omega_power_2",
    "4_Des_F_R_LCDM_2", "5_hdsk_TG_2", "5_hdsk_TG_1", "5_Horava_eta_1", "5_quint_exponential", "4_Des_mc5e_CPL_1",
    "3_RPH_KT_power_1", "5_scg", "5_Kmouflage_3", "2_PEFT_gamma5_exp_2", "2_PEFT_gamma3_lin_2", "4_Des_F_R_wCDM_1",
    "5_ADE_1", "2_PEFT_wDE_wCDM_1",
]
MINI = ["2_PEFT_Omega_power_2", "2_PEFT_wDE_wCDM_1", "3_RPH_KT_power_1", "4_Des_F_R_LCDM_2", "4_Des_mc5e_CPL_1",
        "5_quint_exponential", "5_hdsk_TG_2", "5_scg"]
SHORT = {
    "bto": "EFTCAMB_back_turn_on",
    "pto": "EFTCAMB_turn_on_time",
    "grt": "EFTCAMB_GR_threshold",
    "stt": "EFTCAMB_stability_time",
}


def overrides_for(var: str, ini: Path) -> dict[str, str]:
    if var == "ref":
        return {}
    if var == "AB2":
        return {"accuracy_boost": "2"}
    if var == "mh":
        return {"EFTCAMB_evolve_metric_h": "T"}
    if var == "ub":
        return {"EFTCAMB_use_background": "T"}
    if var == "dphiT":
        return {"EFTCAMB_evolve_delta_phi": "T"}
    if var == "dphiF":
        return {"EFTCAMB_evolve_delta_phi": "F"}
    if var.startswith("ic"):
        return {"EFT_IC_type": var[2:]}
    key = SHORT[var[:3]]
    keys = harness.parse_ini(ini)
    v = float(str(keys.get(key, DEFAULTS[key])).lower().replace("d", "e"))
    fac = 10.0 if var.endswith("x10") else 0.1
    return {key: f"{v * fac:.6e}"}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("variations")
    p.add_argument("--models", default="*.ini")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--tree", default="A0")
    p.add_argument("--outname", default="")
    p.add_argument("--list", default="", help="comma-separated model stems (overrides --models); 'rep' = REPRESENTATIVE")
    a = p.parse_args()
    inis = harness.model_inis(a.models)
    if a.list:
        names = {"rep": REPRESENTATIVE, "core": CORE, "mini": MINI}.get(a.list) or a.list.split(",")
        inis = [Path(n) if "/" in n else harness.PARAMS_DIR / f"{n}.ini" for n in names]
    tasks = []
    # each spec: VAR, or OUTNAME:VAR@TREE (TREE = tree name or path) to run a variation on another tree
    for spec in a.variations.split(","):
        outname, var, tree = (a.outname or spec), spec, a.tree
        if ":" in spec:
            outname, rest = spec.split(":", 1)
            var, tree = rest.split("@", 1) if "@" in rest else (rest, a.tree)
        outdir = WORK / "runs" / outname
        for ini in inis:
            out = outdir / f"{ini.stem}.npz"
            if out.with_suffix(".json").exists():
                continue
            tasks.append((var, ini, out, tree))
    print(f"{len(tasks)} runs", flush=True)

    def one(t):
        var, ini, out, tree = t
        st = harness.run_model(tree, ini, out, overrides=overrides_for(var, ini), threads=2, timeout=1500)
        print(f"{time.strftime('%H:%M:%S')} {var:8s} {ini.stem:32s} {st['status']:16s} {st['elapsed']:7.1f}s",
              flush=True)
        return st

    with ThreadPoolExecutor(a.jobs) as ex:
        list(ex.map(one, tasks))


if __name__ == "__main__":
    main()

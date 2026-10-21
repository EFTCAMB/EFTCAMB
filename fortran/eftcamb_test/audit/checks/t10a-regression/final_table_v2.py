"""Track 10a final per-model table: worst inherited-corrected deviation (double ratio) A0 vs B1 at default settings.

Reads dr_default_v4.json (double_ratio.py runs_A0_itb2 runs_B1_itb1, best matched GR reference per model) and assigns
the classification established by the follow-up experiments (see the track report).  Prints one JSON object per model.
Usage: python final_table_v2.py WORKDIR
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

COLS = [("TT", "TT"), ("EE", "EE"), ("TE", "TE"), ("BB", "BB_lens"), ("PP", "PP"), ("Pk", "Pk")]
UNCONVERGED = {"3_RPH_MB_const_1", "5_Kmimic_2", "5_Horava_SolSyst_1", "5_Horava_SolSyst_2", "5_Horava_SolSyst_3"}


def classify(m: str) -> tuple[str, bool]:
    if m in UNCONVERGED:
        return "known not converged at default accuracy in either tree (accuracy-limited, not classifiable)", False
    if m.startswith(("2_PEFT_Omega", "2_PEFT_gamma3", "2_PEFT_gamma4", "2_PEFT_gamma5", "3_RPH", "4_Des_F_R")):
        return ("inherited/accuracy: converged runs (ab2+lab2+itb2, mnu=0) of Omega_const_1, gamma3/4_const_1, "
                "RPH_MK_const_1, F_R_LCDM_1 give unlensed DR<=2e-4 at ell>30, PP<=2e-4, Pk<=3e-4; residual ell<30 "
                "<=4e-3 equals the codes' own low-ell accuracy; low-ell lensed BB from CAMB lensing method"), True
    if m.startswith("5_Horava_eta"):
        return ("inherited/accuracy: eta_4 converged DR<=3e-3 (TT/EE), its own default-vs-converged error 3e-2; "
                "tBB from CAMB v2 massive-nu tensor change (mnu=0: 1.3e-4)"), True
    if m.startswith("5_hdsk"):
        return ("inherited: C_l interpolation between sampled multipoles near spectral troughs; with "
                "l_sample_boost=50 in both trees TG_1 DR <=1.2e-3 (uEE ell 1713), tBB 1.8e-4 (was 1.3e-2)"), True
    if m.startswith("5_Kmouflage"):
        return ("inherited: multipole interpolation (Kmouflage_2 EE trough 2e-2 -> 2e-3 with l_sample_boost=50) "
                "plus default accuracy (own errors up to 1e-1 at low ell)"), True
    if m.startswith("5_Kmimic"):
        return ("inherited/accuracy: Cl converge (ab2: DR<=2e-3 vs own errors 2-5e-2); P(k) jump at k~0.21/Mpc "
                "unresolved in both trees"), True
    if m == "5_scg":
        return "inherited: massive-neutrino CAMB v2 changes (mnu=0 DR<=5e-4; converged mnu=0 <=4e-4)", True
    return "DR >= 1e-4 unexplained", False


def main() -> None:
    work = Path(sys.argv[1])
    d = json.loads((work / "dr_default_v4.json").read_text())["models"]
    for m, rec in d.items():
        if m == "1_EFT_GR" or "double" not in rec:
            continue
        dd = rec["double"]
        vals = {lab: dd[k]["max"] for lab, k in COLS if k in dd}
        worst = max(vals.values())
        if worst < 1e-4:
            note, ok = "inherited (all double ratios < 1e-4)", True
        else:
            note, ok = classify(m)
        tbb = dd.get("tBB", {}).get("max", 0.0)
        print(json.dumps({"model": m, "gr_ref": rec["gr_ref"], "vals": vals, "tBB": tbb, "worst": worst,
                          "pass": ok, "note": note}))


if __name__ == "__main__":
    main()

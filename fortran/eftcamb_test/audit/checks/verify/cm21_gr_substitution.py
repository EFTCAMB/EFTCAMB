"""Quantify the GR Einstein-equation substitutions hard-coded in the 21cm line-of-sight sources (equations.f90,
output_window_sources, window_21cm branch, line_distortions) when EFTCAMB is active.

The line_distortions source starts with
    s(1) = (sigma/adotoa/3 - etak/adotoa**2/3)*wing*exptau*k ;  s(2) = -dgrho/(6 adotoa**2)*wing*exptau + ...
i.e. k/(3 adotoa) * (sigma - z_GR) * wing * exptau  with  z_GR = (etak + dgrho/(2k))/adotoa  (GR 00 equation).
Kinematically (sigma = (h' + 6 eta')/(2k), z = h'/(2k)) sigma - z = 3*etakdot/k**2 in any theory, so the
theory-independent form is k/(3 adotoa)*(sigma - z) = etakdot/(k adotoa).  The code passes etakdot to the routine but
uses z_GR.  We compare A_code = k(sigma - z_GR)/3 against A_kin = etakdot/k (= k(sigma - z)/3 with the true z) using
the EFTCAMB time-step cache along the evolution (get_scale_evolution), for a GR control and an EFT model.
Also compare the GR sigma' relation sigma'_GR = -2 adotoa sigma + etak (dgpi/k neglected) used in the remaining
terms of s(2) (the 2 etak/Kf terms), against the cached sigmadot.
"""

import json
import os
import sys
from pathlib import Path

import numpy as np

os.environ["OMP_NUM_THREADS"] = "1"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

SCR = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit")
TREE = SCR / "trees" / "A0"
OUT = SCR / "work" / "verify" / "21cm"
OUT.mkdir(parents=True, exist_ok=True)
PARAMS = HERE.parents[2] / "parameters"
KS = [0.005, 0.02, 0.05, 0.1]
MODELS = sys.argv[1:] or ["2_PEFT_wDE_LCDM_1", "2_PEFT_Omega_const_1", "2_PEFT_Omega_lin_1",
                          "4_Des_F_R_LCDM_1", "2_PEFT_gamma1_const_1"]


def analyse(d: dict) -> dict:
    a = d["evol_a"][0] if d["evol_a"].ndim == 2 else d["evol_a"]  # arrays are [k, a]
    kv = d["evol_k"][:, 0] if d["evol_k"].ndim == 2 else d["evol_k"]
    k = kv[:, None]
    f = {n: d[f"evol_{n}"] for n in ("etak", "etakdot", "z", "sigma", "sigmadot", "dgrho", "adotoa")}
    zgr = (f["etak"] + 0.5 * f["dgrho"] / k) / f["adotoa"]
    a_code = k * (f["sigma"] - zgr) / 3
    a_kin = f["etakdot"] / k
    a_true = k * (f["sigma"] - f["z"]) / 3
    sdot_gr = -2 * f["adotoa"] * f["sigma"] + f["etak"]
    res = {}
    for a0 in (1 / 31.0, 0.1, 0.5, 1.0):
        i = int(np.argmin(abs(a - a0)))
        row = {}
        for j, kk in enumerate(kv):
            row[f"k={kk:g}"] = {
                "z_true": float(f["z"][j, i]), "z_GR": float(zgr[j, i]),
                "rel_z": float(abs(zgr[j, i] - f["z"][j, i]) / abs(f["z"][j, i])),
                "A_code/A_kin": float(a_code[j, i] / a_kin[j, i]),
                "A_true/A_kin": float(a_true[j, i] / a_kin[j, i]),
                "sigmadot_GR/sigmadot": float(sdot_gr[j, i] / f["sigmadot"][j, i]),
            }
        res[f"a={a[i]:.4g}"] = row
    return res


if __name__ == "__main__":
    summary = {}
    for model in MODELS:
        npz = OUT / f"evol_{model}.npz"
        if not npz.exists():
            st = harness.run_model(str(TREE), str(PARAMS / f"{model}.ini"), str(npz), evol_k=KS, threads=1,
                                   timeout=900)
            if st["status"] != "ok":
                print(model, st["status"], st.get("log_tail"))
                continue
        d = harness.load(npz)
        if "evol_etak" not in d:
            print(model, "no evolution output", d["meta"].get("errors") if isinstance(d["meta"], dict) else d["meta"])
            continue
        summary[model] = analyse(d)
        print(f"== {model}")
        for ak, row in summary[model].items():
            for kk, v in row.items():
                print(f"  {ak:9s} {kk:8s} rel|z_GR-z|/|z|={v['rel_z']:.3e}  A_code/A_kin={v['A_code/A_kin']:+.4f}  "
                      f"A_true/A_kin={v['A_true/A_kin']:+.6f}  sigmadot_GR/sigmadot={v['sigmadot_GR/sigmadot']:+.4f}")
    (OUT / "cm21_gr_substitution.json").write_text(json.dumps(summary, indent=1))

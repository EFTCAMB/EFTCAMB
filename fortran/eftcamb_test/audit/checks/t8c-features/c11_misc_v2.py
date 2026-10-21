"""Check 11: miscellaneous new CAMB v2 features with EFT on (designer f(R) LCDM_1 vs GR):
  heat     : Reion.include_heating=True (new reionization heating switch) -> relative P(k) effect EFT vs GR
  accBB    : Accuracy.AccurateBB=True (new full Gauss-Legendre curved-sky lensing method) -> runs, BB change vs default
  spk      : SPkNonLinear wrapping Halofit (new non-linear model) -> runs, finite
  margin   : lens_output_margin=400 (new first-class parameter) -> lensed Cls to lmax unchanged at the 1e-3 level
"""

from __future__ import annotations

import json
import sys

import numpy as np
from t8c_lib import WORK, cl_diffs, load, model_ini, pk_diff, run_many

OUT = WORK / "c11"
OV = {"get_tensor_cls": "F", "feedback_level": "0", "num_redshiftwindows": "0", "transfer_kmax": "10"}
BASE = {"set_for_lmax": {"lmax": 2000, "lens_potential_accuracy": 1, "nonlinear": False},
        "matter_power": {"redshifts": [0.0], "kmax": 10.0}, "outputs": ["cmb", "pk_lin"], "cl_lmax": 2000}
VARIANTS = {
    "base": {},
    "heat": {"set_attrs": {"Reion.include_heating": True}},
    "accBB": {"set_attrs": {"Accuracy.AccurateBB": True}},
    "margin": {"set_for_lmax": {"lmax": 2000, "lens_potential_accuracy": 1, "nonlinear": False,
                                "lens_output_margin": 400}},
}
MODELS = {"GR": "1_EFT_GR", "fR": "4_Des_F_R_LCDM_1", "g3": "2_PEFT_gamma3_const_1"}


def jobs():
    js = []
    for tag, m in MODELS.items():
        for v, a in VARIANTS.items():
            js.append(dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{tag}_{v}.npz", overrides=OV,
                           args=BASE | a))
        js.append(dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{tag}_spk.npz", overrides=OV,
                       args=BASE | {"set_classes_nl": "SPkNonLinear", "NonLinear": "NonLinear_pk",
                                    "outputs": ["pk_lin", "pk_nl"]}))
    return js


def report():
    res = {}
    for tag in MODELS:
        try:
            b = load(OUT / f"{tag}_base.npz")
        except FileNotFoundError:
            res[tag] = "missing"
            continue
        row = {}
        for v in ("heat", "accBB", "margin"):
            try:
                a = load(OUT / f"{tag}_{v}.npz")
            except FileNotFoundError:
                st = json.loads((OUT / f"{tag}_{v}.json").read_text())
                row[v] = st["status"] + " " + str(st.get("log_tail", [])[-2:])
                continue
            row[v] = {"cl": cl_diffs(a["cl_total"], b["cl_total"], lmax=2000), "pk": pk_diff(a, b, kmax=10)}
            if v == "heat":
                kk = a["pk_k"]
                r = a["pk"][0] / b["pk"][0] - 1
                row[v]["pk_rel_at_k"] = {f"{k:g}": float(np.interp(k, kk, r)) for k in (0.1, 1.0, 5.0)}
        try:
            s = load(OUT / f"{tag}_spk.npz")
            rr = s["pknl"] / s["pk"]
            row["spk"] = {"finite": bool(np.all(np.isfinite(rr))), "ratio_min": float(rr.min()),
                          "ratio_max": float(rr.max())}
        except FileNotFoundError:
            st = json.loads((OUT / f"{tag}_spk.json").read_text())
            row["spk"] = st["status"] + " " + str(st.get("log_tail", [])[-2:])
        res[tag] = row
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        for r in run_many(jobs(), nproc=int(sys.argv[2]) if len(sys.argv) > 2 else 4):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report()

"""Spectra impact of the Horava Hdot fix, with an accuracy-noise control (verify track).

Runs (harness, threads=1): patched tree P_horava (t4b-derivs, differs from A0 only in 008p1_Horava.f90) at default
accuracy for the Horava test models, and A0/P_horava at accuracy_boost=2,l_accuracy_boost=2 for eta_4/SolSyst_3/eta_2.
Compares A0(default ref) vs P(default), A0(default) vs A0(boost) [noise], A0(boost) vs P(boost).

Usage: python horava_hdot_spectra.py run|compare
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
P = AUDIT / "work" / "t4b-derivs" / "P_horava"
A0 = AUDIT / "trees" / "A0"
OUT = AUDIT / "work" / "verify" / "horava_spectra"
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
MODELS = ["5_Horava_eta_1", "5_Horava_eta_2", "5_Horava_eta_3", "5_Horava_eta_4", "5_Horava_SolSyst_3"]
BOOST_MODELS = ["5_Horava_eta_4", "5_Horava_SolSyst_3", "5_Horava_eta_2"]
BOOST = {"accuracy_boost": "2", "l_accuracy_boost": "2"}


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    inis = [PARAMS / f"{m}.ini" for m in MODELS]
    st = harness.run_batch(P, inis, OUT / "P_def", jobs=4, threads=1)
    print(json.dumps([(s["model_ini"].split("/")[-1], s["status"], s.get("elapsed")) for s in st]))
    for tree, tag in ((A0, "A0_boost"), (P, "P_boost")):
        st = harness.run_batch(tree, [PARAMS / f"{m}.ini" for m in BOOST_MODELS], OUT / tag, jobs=3, threads=1,
                               overrides=BOOST)
        print(tag, json.dumps([(s["model_ini"].split("/")[-1], s["status"], s.get("elapsed")) for s in st]))


def lowl_tt(a, b, lmax=30):
    ta, tb = a["cl_total"][2:lmax + 1, 0], b["cl_total"][2:lmax + 1, 0]
    return float(np.max(np.abs(ta - tb) / np.abs(ta)))


def summary(a, b):
    c = harness.compare(a, b)
    da, db = harness.load(a), harness.load(b)
    tt_hi = float(np.max(np.abs(da["cl_total"][200:2501, 0] - db["cl_total"][200:2501, 0]) /
                         np.abs(da["cl_total"][200:2501, 0])))
    return {"TT_l2_30": lowl_tt(da, db), "TT_l200_2500": tt_hi, "EE": c["cl"]["cl_total"]["EE"] if "cl_total" in
            c["cl"] else None, "worst_cl": c["worst_cl"], "worst_pk": c["worst_pk"],
            "sigma8": [float(np.ravel(da["sigma8"])[-1]), float(np.ravel(db["sigma8"])[-1])]}


def compare():
    res = {}
    for m in MODELS:
        res[f"{m}: A0def vs Pdef"] = summary(AUDIT / "runs" / "A0" / f"{m}.npz", OUT / "P_def" / f"{m}.npz")
    for m in BOOST_MODELS:
        res[f"{m}: A0def vs A0boost (noise)"] = summary(AUDIT / "runs" / "A0" / f"{m}.npz", OUT / "A0_boost" / f"{m}.npz")
        res[f"{m}: Pdef vs Pboost (noise)"] = summary(OUT / "P_def" / f"{m}.npz", OUT / "P_boost" / f"{m}.npz")
        res[f"{m}: A0boost vs Pboost"] = summary(OUT / "A0_boost" / f"{m}.npz", OUT / "P_boost" / f"{m}.npz")
    for k, v in res.items():
        print(f"{k:45s} TT(2-30) {v['TT_l2_30']:.3e}  TT(200-2500) {v['TT_l200_2500']:.3e}  worst_cl {v['worst_cl']}"
              f"  worst_pk {v['worst_pk']}  sigma8 {v['sigma8'][0]:.5f}->{v['sigma8'][1]:.5f}")
    (OUT / "summary.json").write_text(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    {"run": run, "compare": compare}[sys.argv[1]]()

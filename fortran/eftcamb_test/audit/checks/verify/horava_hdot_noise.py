"""Accuracy-noise control for finding horava-hdot-inconsistent (verify track).

Runs boosted-accuracy versions of selected Horava models on A0 and on the patched tree P_horava (t4b-derivs; differs
from A0 only in 008p1_Horava.f90 compute_H_derivs) and compares:
  noise  = A0(default) vs A0(boost)       [what the accuracy setting alone changes]
  effect = A0(boost)   vs P(boost)        [what the Hdot fix changes at higher accuracy]

Usage: python horava_hdot_noise.py run TAG MODEL OVERRIDES_JSON [threads]   (TAG chooses the output subdirectory)
       python horava_hdot_noise.py compare
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

AUDIT = harness.AUDIT_ROOT
P = AUDIT / "work" / "t4b-derivs" / "P_horava"
A0 = AUDIT / "trees" / "A0"
OUT = AUDIT / "work" / "verify" / "horava_spectra"
PARAMS = harness.PARAMS_DIR


def run(jobs: list[tuple[str, str, str, dict, int]]):
    """jobs: (tree_tag, subdir, model, overrides, threads)."""
    trees = {"A0": A0, "P": P}

    def one(j):
        tag, sub, model, ov, th = j
        (OUT / sub).mkdir(parents=True, exist_ok=True)
        return harness.run_model(trees[tag], PARAMS / f"{model}.ini", OUT / sub / f"{model}.npz", overrides=ov,
                                 threads=th, timeout=560)

    with ThreadPoolExecutor(len(jobs)) as ex:
        for j, s in zip(jobs, ex.map(one, jobs)):
            print(j[:3], s["status"], s.get("elapsed"))


def tt_rel(da, db, l0, l1):
    ta, tb = da["cl_total"][l0:l1 + 1, 0], db["cl_total"][l0:l1 + 1, 0]
    return float(np.max(np.abs(ta - tb) / np.abs(ta)))


def ee_rel(da, db, l0, l1):
    ta, tb = da["cl_total"][l0:l1 + 1, 1], db["cl_total"][l0:l1 + 1, 1]
    return float(np.max(np.abs(ta - tb) / np.abs(ta)))


def summary(a, b):
    da, db = harness.load(a), harness.load(b)
    return {"TT_2_30": tt_rel(da, db, 2, 30), "TT_200_2500": tt_rel(da, db, 200, 2500),
            "EE_30_2000": ee_rel(da, db, 30, 2000),
            "sigma8": [float(np.ravel(da["sigma8"])[-1]), float(np.ravel(db["sigma8"])[-1])]}


def compare():
    pairs = []
    for m in ["5_Horava_eta_2"]:
        pairs += [(f"{m} noise A0def vs A0boost2", AUDIT / "runs" / "A0" / f"{m}.npz", OUT / "A0_boost" / f"{m}.npz"),
                  (f"{m} effect A0boost2 vs Pboost2", OUT / "A0_boost" / f"{m}.npz", OUT / "P_boost" / f"{m}.npz"),
                  (f"{m} effect A0def vs Pdef", AUDIT / "runs" / "A0" / f"{m}.npz", OUT / "P_def" / f"{m}.npz")]
    for m in ["5_Horava_eta_4", "5_Horava_SolSyst_3"]:
        pairs += [(f"{m} noise A0def vs A0b15", AUDIT / "runs" / "A0" / f"{m}.npz", OUT / "A0_b15" / f"{m}.npz"),
                  (f"{m} noise Pdef vs Pb15", OUT / "P_def" / f"{m}.npz", OUT / "P_b15" / f"{m}.npz"),
                  (f"{m} effect A0b15 vs Pb15", OUT / "A0_b15" / f"{m}.npz", OUT / "P_b15" / f"{m}.npz"),
                  (f"{m} effect A0def vs Pdef", AUDIT / "runs" / "A0" / f"{m}.npz", OUT / "P_def" / f"{m}.npz")]
    res = {}
    for name, a, b in pairs:
        if not (Path(a).exists() and Path(b).exists()):
            print(f"{name:48s} missing")
            continue
        v = summary(a, b)
        res[name] = v
        print(f"{name:48s} TT(2-30) {v['TT_2_30']:.3e}  TT(200-2500) {v['TT_200_2500']:.3e}  "
              f"EE(30-2000) {v['EE_30_2000']:.3e}  sigma8 {v['sigma8'][0]:.5f}->{v['sigma8'][1]:.5f}")
    (OUT / "noise_summary.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "compare":
        compare()
    else:
        run(json.loads(sys.argv[2]))

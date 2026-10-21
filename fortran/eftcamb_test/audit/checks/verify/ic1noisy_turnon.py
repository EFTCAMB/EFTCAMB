"""Verifier (ic1-pidot-transient-noisy-coefficients): end-to-end attribution runs for 5_Kmimic_1 (and TG) on A0.

Questions:
  (a) Is the K-mimic turn-on dependence tied to where a_on sits relative to the nodes of the model's background
      interpolation table (default model_background_num_points = 1000 on ln a in [ln 1e-9, 0], so 1e-8 and 1e-7 are
      exactly nodes)?  -> tiny shifts of a_on and different num_points.
  (b) Is the turn-on 1e-9 blow-up a table-edge effect (a_on = model_background_a_ini = 1e-9)? -> a_ini = 1e-11.
  (c) (EFT_IC_type = 0 is rejected by the code: "Invalid EFT_IC_type"; dropped.)
Usage: python ic1noisy_turnon.py [group ...]   groups: kmimic, kmimic2, tg
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy" / "turnon"
PIDFIX = harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_pidfix"
KS = (0.1, 0.5, 1.0, 2.0)

KM = "5_Kmimic_1"
GROUPS = {
    "kmimic": [
        ("GR", "A0", "1_EFT_GR", {}),
        ("km_ic1_on1e-8", "A0", KM, {"EFTCAMB_turn_on_time": "1e-8"}),
        ("km_ic1_on1.003e-8", "A0", KM, {"EFTCAMB_turn_on_time": "1.003e-8"}),
        ("km_ic1_on1.01e-8", "A0", KM, {"EFTCAMB_turn_on_time": "1.01e-8"}),
        ("km_ic1_on1.03e-8", "A0", KM, {"EFTCAMB_turn_on_time": "1.03e-8"}),
        ("km_ic1_on0.97e-8", "A0", KM, {"EFTCAMB_turn_on_time": "0.97e-8"}),
        ("km_ic1_on3e-8", "A0", KM, {"EFTCAMB_turn_on_time": "3e-8"}),
        ("km_ic1_on1e-7", "A0", KM, {"EFTCAMB_turn_on_time": "1e-7"}),
        ("km_ic1_on1.01e-7", "A0", KM, {"EFTCAMB_turn_on_time": "1.01e-7"}),
        ("km_ic1_on1e-8_np1500", "A0", KM, {"EFTCAMB_turn_on_time": "1e-8", "model_background_num_points": "1500"}),
        ("km_ic1_on1e-8_np5000", "A0", KM, {"EFTCAMB_turn_on_time": "1e-8", "model_background_num_points": "5000"}),
        ("km_ic1_on1e-9", "A0", KM, {"EFTCAMB_turn_on_time": "1e-9"}),
        ("km_ic1_on1e-9_ai1e-11", "A0", KM, {"EFTCAMB_turn_on_time": "1e-9", "model_background_a_ini": "1e-11"}),
        ("km_ic1_on1e-10_ai1e-11", "A0", KM, {"EFTCAMB_turn_on_time": "1e-10", "model_background_a_ini": "1e-11"}),
    ],
    "kmimic2": [
        ("km_pidfix_on1e-8", str(PIDFIX), KM, {"EFTCAMB_turn_on_time": "1e-8"}),
        ("km_pidfix_on1.01e-8", str(PIDFIX), KM, {"EFTCAMB_turn_on_time": "1.01e-8"}),
        ("km_pidfix_on1e-8_np5000", str(PIDFIX), KM, {"EFTCAMB_turn_on_time": "1e-8",
                                                      "model_background_num_points": "5000"}),
        ("km_B1_on1e-8", "B1", KM, {"EFTCAMB_turn_on_time": "1e-8"}),
        ("km_B1_on1e-7", "B1", KM, {"EFTCAMB_turn_on_time": "1e-7"}),
        ("km_B1_on1.01e-8", "B1", KM, {"EFTCAMB_turn_on_time": "1.01e-8"}),
        ("GR_B1", "B1", "1_EFT_GR", {}),
    ],
    "tg": [
        ("tg1_A0", "A0", "5_hdsk_TG_1", {}),
        ("tg2_A0", "A0", "5_hdsk_TG_2", {}),
        ("tg2_pidfix", str(PIDFIX), "5_hdsk_TG_2", {}),
        ("tg2_A0_on1.01e-8", "A0", "5_hdsk_TG_2", {"EFTCAMB_turn_on_time": "1.01e-8"}),
        ("tg1_A0_on1.01e-8", "A0", "5_hdsk_TG_1", {"EFTCAMB_turn_on_time": "1.01e-8"}),
        ("tg2_A0_on2e-8", "A0", "5_hdsk_TG_2", {"EFTCAMB_turn_on_time": "2e-8"}),
        ("tg1_A0_on2e-8", "A0", "5_hdsk_TG_1", {"EFTCAMB_turn_on_time": "2e-8"}),
        ("tg1_B1", "B1", "5_hdsk_TG_1", {}),
        ("tg2_B1", "B1", "5_hdsk_TG_2", {}),
        ("tg2_B1_on1.01e-8", "B1", "5_hdsk_TG_2", {"EFTCAMB_turn_on_time": "1.01e-8"}),
        ("tg1_B1_on1.01e-8", "B1", "5_hdsk_TG_1", {"EFTCAMB_turn_on_time": "1.01e-8"}),
    ],
}


def job(args):
    name, tree, model, ov = args
    out = OUT / f"{name}.npz"
    if not out.exists():
        st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=ov, threads=1, timeout=900)
        if st["status"] != "ok":
            return name, st["status"]
    return name, harness.load(out)


def pk_at(d, k):
    return float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))


def cl_diff(a, b):
    out = []
    for j in (0, 1):
        x, y = a["cl_total"][2:2501, j], b["cl_total"][2:2501, j]
        r = np.abs(x / y - 1)
        out.append((float(r.max()), int(np.argmax(r)) + 2))
    return out


def main():
    groups = sys.argv[1:] or ["kmimic"]
    jobs = [j for g in groups for j in GROUPS[g]]
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(job, jobs))
    for name, d in res.items():
        if isinstance(d, str):
            print(f"{name}: {d}")
    if "kmimic" in groups or "kmimic2" in groups:
        grs = {"A0": res.get("GR") or harness.load(OUT / "GR.npz")}
        if "GR_B1" in res:
            grs["B1"] = res["GR_B1"]
        print("5_Kmimic_1 P(k)/P_GR - 1 at k = " + ", ".join(f"{k:g}" for k in KS) + " /Mpc, sigma8")
        for name, d in res.items():
            if not name.startswith("km") or isinstance(d, str):
                continue
            gr = grs["B1"] if "_B1_" in name and "B1" in grs else grs["A0"]
            vals = "  ".join(f"{pk_at(d, k) / pk_at(gr, k) - 1:+.4e}" for k in KS)
            print(f"  {name:28s} {vals}   s8 {float(d['sigma8'][0]):.5f}")
    if "tg" in groups:
        pairs = [("tg2_A0", "tg1_A0"), ("tg2_pidfix", "tg1_A0"),
                 ("tg2_A0_on1.01e-8", "tg1_A0_on1.01e-8"), ("tg2_A0_on2e-8", "tg1_A0_on2e-8"),
                 ("tg2_B1", "tg1_B1"), ("tg2_B1_on1.01e-8", "tg1_B1_on1.01e-8"), ("tg1_A0_on1.01e-8", "tg1_A0"),
                 ("tg1_A0_on2e-8", "tg1_A0")]
        for p, q in pairs:
            a, b = res.get(p), res.get(q)
            if a is None or b is None or isinstance(a, str) or isinstance(b, str):
                continue
            (tt, ltt), (ee, lee) = cl_diff(a, b)
            print(f"  {p:20s} vs {q:20s}: max|dTT/TT| = {tt:.2e} (ell {ltt}), max|dEE/EE| = {ee:.2e} (ell {lee})")


if __name__ == "__main__":
    main()

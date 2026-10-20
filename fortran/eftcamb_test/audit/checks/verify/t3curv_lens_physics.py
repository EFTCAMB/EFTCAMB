"""Adversarial physics-lens verifier for finding t3-eft-curvature-ignored.

All runs via harness.run_model (threads=1, <=4 concurrent). Outputs in work/verify/t3curv_lens/.
Reference is always plain upstream CAMB v2 (B0, EFTflag=0) at the same cosmology.

Checks
  S  scaling of the GR-limit discrepancy with omk (Pure EFT, Omega=0, w=-1), including Planck-sized |omk|~1e-3
  W  second GR-limit pure-EFT model (EFTwDE=1, w0=-1) at omk=0.05
  R  designer RPH with all alpha functions zero (GR limit) at omk=0 and 0.05: spectra AND background H(z)
  T  tensors: GR-limit pure EFT vs B0 at omk=0.05 (derivst also takes adotoa from compute_adotoa)
  D  attribution: diagnostic tree (copy of A0 where derivs adds grhok to adotoa, uses CAMB's adotdota and the curved
     cothxor when not flat) must bring the GR-limit model back to B0 if the stated root cause is the whole story
  L  ell dependence of the omk=0.05 TT ratio
Usage: python t3curv_lens_physics.py [run|analyse]
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "t3curv_lens"
DIAG = harness.AUDIT_ROOT / "work" / "verify" / "curvdiag_tree"
P = harness.PARAMS_DIR
PEFT = P / "2_PEFT_Omega_const_1.ini"
PEFTW = P / "2_PEFT_wDE_LCDM_1.ini"
RPH = P / "3_RPH_MK_const_1.ini"
GR = P / "1_EFT_GR.ini"
COMMON = {"num_redshiftwindows": "0", "get_tensor_cls": "F", "l_max_scalar": "2200"}
TENS = {"get_tensor_cls": "T", "l_max_tensor": "1500", "k_eta_max_tensor": "3000"}
OMKS_S = ["0", "0.001", "-0.001", "0.003", "0.01", "-0.01", "0.05", "-0.05"]


def jobs():
    j = []
    for omk in OMKS_S:
        ov = {**COMMON, "omk": omk}
        j.append(("A0", PEFT, {**ov, "PureEFTmodelOmega": "0"}, f"A0_PEFT0_omk{omk}"))
        j.append(("B0", GR, ov, f"B0_GR_omk{omk}"))
    for omk in ("0", "0.05"):
        ov = {**COMMON, "omk": omk}
        j.append(("A0", PEFTW, ov, f"A0_PEFTW_omk{omk}"))
        j.append(("A0", RPH, {**ov, "RPHmassPmodel": "0", "RPHkineticitymodel": "0"}, f"A0_RPH0_omk{omk}"))
        j.append(("A0", PEFT, {**ov, **TENS, "PureEFTmodelOmega": "0"}, f"A0_PEFT0T_omk{omk}"))
        j.append(("B0", GR, {**ov, **TENS}, f"B0_GRT_omk{omk}"))
    for omk in ("0", "0.01", "0.05", "-0.05"):
        ov = {**COMMON, "omk": omk}
        j.append((str(DIAG), PEFT, {**ov, "PureEFTmodelOmega": "0"}, f"D_PEFT0_omk{omk}"))
    j.append((str(DIAG), GR, {**COMMON, "omk": "0.05"}, "D_GR_omk0.05"))
    return j


def run_all():
    OUT.mkdir(parents=True, exist_ok=True)
    todo = [x for x in jobs() if not (OUT / f"{x[3]}.npz").exists()]

    def one(x):
        tree, ini, ov, name = x
        return harness.run_model(tree, ini, OUT / f"{name}.npz", overrides=ov, threads=1, timeout=900)

    with ThreadPoolExecutor(4) as ex:
        for s in ex.map(one, todo):
            if s["status"] != "ok":
                print("FAILED", s["out"], s["status"], s.get("log_tail", [])[-4:])


def row(a, b):
    r = harness.compare(OUT / f"{a}.npz", OUT / f"{b}.npz")
    t = r["cl"].get("total", {})
    A, B = harness.load(OUT / f"{a}.npz"), harness.load(OUT / f"{b}.npz")
    tt_a, tt_b = A["cl_total"][2:2201, 0], B["cl_total"][2:2201, 0]
    rel = np.abs(tt_a / tt_b - 1)
    return {
        "TT": t.get("TT"), "TT_l>=30": float(rel[28:].max()), "EE": t.get("EE"), "TE": t.get("TE"),
        "BB": t.get("BB"), "PP": r["cl"].get("lens_potential", {}).get("PP"), "tensorTT": r["cl"].get("tensor", {}).get("TT"),
        "tensorBB": r["cl"].get("tensor", {}).get("BB"), "pk": r["worst_pk"], "sigma8": r.get("sigma8"),
        "bg_H": r.get("bg_H"), "TT_l2_ratio": float(tt_a[0] / tt_b[0]),
    }


def analyse():
    res = {}
    for omk in OMKS_S:
        res[f"S PEFT0 vs B0 omk={omk}"] = row(f"A0_PEFT0_omk{omk}", f"B0_GR_omk{omk}")
    for omk in ("0", "0.05"):
        res[f"W PEFT(EFTwDE=1,w0=-1) vs B0 omk={omk}"] = row(f"A0_PEFTW_omk{omk}", f"B0_GR_omk{omk}")
        res[f"R RPH-zero vs B0 omk={omk}"] = row(f"A0_RPH0_omk{omk}", f"B0_GR_omk{omk}")
        res[f"T PEFT0+tensors vs B0 omk={omk}"] = row(f"A0_PEFT0T_omk{omk}", f"B0_GRT_omk{omk}")
    for omk in ("0", "0.01", "0.05", "-0.05"):
        res[f"D diag-patched PEFT0 vs B0 omk={omk}"] = row(f"D_PEFT0_omk{omk}", f"B0_GR_omk{omk}")
    res["D diag-patched GR vs B0 omk=0.05 (inert for EFTflag=0)"] = row("D_GR_omk0.05", "B0_GR_omk0.05")
    A, B = harness.load(OUT / "A0_PEFT0_omk0.05.npz"), harness.load(OUT / "B0_GR_omk0.05.npz")
    res["L TT ratio A0/B0 omk=0.05 by ell"] = {
        str(ell): float(A["cl_total"][ell, 0] / B["cl_total"][ell, 0]) for ell in (2, 3, 4, 5, 8, 10, 20, 30, 50, 100,
                                                                                  300, 1000, 2000)
    }
    R = harness.load(OUT / "A0_RPH0_omk0.05.npz")
    Bz = harness.load(OUT / "B0_GR_omk0.05.npz")
    res["R RPH-zero H(z)/H_B0(z)-1 at z grid"] = {
        f"{z:.3g}": float(h / hb - 1) for z, h, hb in zip(R["bg_z"], R["bg_H"], Bz["bg_H"])
    }
    (OUT / "summary.json").write_text(json.dumps(res, indent=1))
    keys = ("TT", "TT_l>=30", "EE", "TE", "PP", "pk", "sigma8", "bg_H", "TT_l2_ratio", "tensorTT", "tensorBB")
    print(f"{'case':58s} " + " ".join(f"{k[:9]:>9s}" for k in keys))
    for lab, r in res.items():
        if lab.startswith(("L ", "R RPH-zero H(z)")):
            print(lab, json.dumps({k: round(v, 5) for k, v in r.items()}))
            continue
        print(f"{lab:58s} " + " ".join(f"{r[k]:9.2e}" if r.get(k) is not None else f"{'-':>9s}" for k in keys))


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == "run":
        run_all()
    analyse()

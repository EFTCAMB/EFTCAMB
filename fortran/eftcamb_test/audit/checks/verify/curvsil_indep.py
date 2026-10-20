"""Independent check of finding curvature-silently-wrong-with-eft (t8c-features).

Claim: omk /= 0 with EFTflag /= 0 is silently accepted and gives wrong spectra even for an exact-LCDM EFT model
(2_PEFT_wDE_LCDM_1), because the perturbation background in derivs ignores curvature.

This script (different method from the repro, which calls camb directly):
  runs     : harness runs (subprocesses, one tree each) at omk in {0, +0.02, -0.02} on
             B0 (upstream CAMB v2, 1_EFT_GR ini = plain LCDM), A0 EFTflag=0, A0 PEFT-LCDM, B1 PEFT-LCDM / B1 GR,
             so the "which side is wrong" question is answered against upstream CAMB, not against EFTCAMB itself.
             The A0 PEFT runs also store get_scale_evolution (derivs path) for two k.
  analyze  : TT/EE/TE/pk/sigma8 ratios and the background check: conformal Hubble rate used in derivs
             (evol_adotoa, EFTCAMB_printing path through derivs) vs the true one a*H(z) from CAMB's background.
Usage: python curvsil_indep.py runs|analyze
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "curvsil"
P = harness.PARAMS_DIR
GR = P / "1_EFT_GR.ini"
PEFT = P / "2_PEFT_wDE_LCDM_1.ini"
COMMON = {"num_redshiftwindows": "0", "get_tensor_cls": "F", "l_max_scalar": "2200"}
OMKS = ("0", "0.02", "-0.02")
EVOL_K = [2e-4, 1e-3]
C_KMS = 299792.458


def jobs():
    j = []
    for omk in OMKS:
        ov = dict(COMMON, omk=omk)
        j.append(("B0", GR, f"B0_GR_omk{omk}", ov, None))
        j.append(("A0", GR, f"A0_GR_omk{omk}", ov, None))
        j.append(("A0", PEFT, f"A0_PEFT_omk{omk}", ov, EVOL_K))
    for omk in ("0", "0.02"):
        ov = dict(COMMON, omk=omk)
        j.append(("B1", GR, f"B1_GR_omk{omk}", ov, None))
        j.append(("B1", PEFT, f"B1_PEFT_omk{omk}", ov, None))
    return j


def runs():
    OUT.mkdir(parents=True, exist_ok=True)

    def one(job):
        tree, ini, name, ov, ek = job
        if (OUT / f"{name}.npz").exists():
            return name, "cached"
        r = harness.run_model(tree, ini, OUT / f"{name}.npz", overrides=ov, evol_k=ek, threads=1, timeout=900)
        return name, r["status"]

    with ThreadPoolExecutor(4) as ex:
        for name, st in ex.map(one, jobs()):
            print(name, st, flush=True)


def tt(d, key="cl_total"):
    return d[key]


def analyze():
    L = {n: harness.load(OUT / f"{n}.npz") for n in [j[2] for j in jobs()] if (OUT / f"{n}.npz").exists()}
    rows = {}

    def ratio(a, b, label):
        A, B = L[a], L[b]
        ca, cb = A["cl_total"], B["cl_total"]
        n = min(len(ca), len(cb), 2201)
        r = ca[2:n] / cb[2:n] - 1
        te = (ca[2:n, 3] - cb[2:n, 3]) / np.sqrt(cb[2:n, 0] * cb[2:n, 1])
        rec = {
            "TT_l2": r[0, 0], "TT_l10": r[8, 0], "TT_l30": r[28, 0], "TT_max_l>=30": np.max(np.abs(r[28:, 0])),
            "EE_max": np.max(np.abs(r[:, 1])), "TE_max_norm": np.max(np.abs(te)),
            "sigma8": float(A["sigma8"][-1] / B["sigma8"][-1] - 1),
            "pk_z0_max": float(np.max(np.abs(A["pk"][-1] / B["pk"][-1] - 1))) if A["pk"].shape == B["pk"].shape else None,
            "bgH_max": float(np.max(np.abs(A["bg_H"] / B["bg_H"] - 1))),
        }
        rows[label] = {k: (float(v) if v is not None else None) for k, v in rec.items()}

    for omk in OMKS:
        ratio(f"A0_GR_omk{omk}", f"B0_GR_omk{omk}", f"A0 EFTflag=0 vs B0 CAMB, omk={omk}")
        ratio(f"A0_PEFT_omk{omk}", f"B0_GR_omk{omk}", f"A0 PEFT-LCDM vs B0 CAMB, omk={omk}")
        ratio(f"A0_PEFT_omk{omk}", f"A0_GR_omk{omk}", f"A0 PEFT-LCDM vs A0 EFTflag=0, omk={omk}")
    for omk in ("0", "0.02"):
        if f"B1_PEFT_omk{omk}" in L and f"B1_GR_omk{omk}" in L:
            ratio(f"B1_PEFT_omk{omk}", f"B1_GR_omk{omk}", f"B1 PEFT-LCDM vs B1 EFTflag=0, omk={omk}")
    hdr = f"{'case':48s} {'TT_l2':>10s} {'TT_l10':>10s} {'TT_l30':>10s} {'TTl>=30':>9s} {'EE_max':>9s} " \
          f"{'TEnorm':>9s} {'sigma8':>10s} {'pk_z0':>9s} {'bgH':>9s}"
    print(hdr)
    for k, v in rows.items():
        print(f"{k:48s} {v['TT_l2']:+10.3e} {v['TT_l10']:+10.3e} {v['TT_l30']:+10.3e} {v['TT_max_l>=30']:9.2e} "
              f"{v['EE_max']:9.2e} {v['TE_max_norm']:9.2e} {v['sigma8']:+10.3e} {v['pk_z0_max'] or 0:9.2e} "
              f"{v['bgH_max']:9.2e}")

    # background check: H used inside derivs vs CAMB's background
    bg = {}
    for omk in OMKS:
        d = L[f"A0_PEFT_omk{omk}"]
        a = d["evol_a"][0]  # the cache field "a" overwrites the harness' evol_a; arrays are [k, time]
        z = 1 / a - 1
        zb, Hb = d["bg_z"], d["bg_H"]  # km/s/Mpc
        sel = (z >= zb[0]) & (z <= zb[-1])
        H_true = a[sel] * np.interp(z[sel], zb, Hb) / C_KMS  # conformal H in 1/Mpc
        H_derivs = d["evol_adotoa"][0, sel]
        # flat-universe prediction of the deficit: adotoa_flat^2 = H_true^2 - grhok/3 = H_true^2 + K, K=-omk H0^2
        H0 = Hb[0] / C_KMS
        H_flat_pred = np.sqrt(H_true**2 - float(omk) * H0**2)
        H_eftfunc = np.interp(a[sel], d["eft_a_grid"], d["eft_adotoa"])
        bg[omk] = {
            "max|adotoa_derivs/H_true-1|": float(np.max(np.abs(H_derivs / H_true - 1))),
            "adotoa_derivs/H_true-1 at a=1": float(H_derivs[-1] / H_true[-1] - 1),
            "max|adotoa_derivs/H_flat_pred-1|": float(np.max(np.abs(H_derivs / H_flat_pred - 1))),
            "max|eft_functions adotoa/H_true-1|": float(np.max(np.abs(H_eftfunc / H_true - 1))),
            "a range": [float(a[sel][0]), float(a[sel][-1])],
        }
    print("\nconformal Hubble rate in derivs (get_scale_evolution) vs CAMB background a*H(z):")
    for k, v in bg.items():
        print(f"  omk={k:6s}", json.dumps(v))
    (OUT / "curvsil_indep_summary.json").write_text(json.dumps({"spectra": rows, "background": bg}, indent=1))


if __name__ == "__main__":
    {"runs": runs, "analyze": analyze}[sys.argv[1]]()

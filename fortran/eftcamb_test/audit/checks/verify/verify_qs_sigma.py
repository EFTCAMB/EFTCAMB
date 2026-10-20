"""Adversarial verification of finding qs-sigma-copy-of-qs-mu (QS_sigma is a verbatim copy of QS_mu).

Independent methods (different code paths from the repro, which used camb.read_ini + get_scale_evolution):
  B) direct ctypes calls of the three exported Fortran symbols
     __eftcamb_abstract_model_MOD_eftcambmodelcomputeqs{mu,sigma,musigma} in A0/camb/camblib.so on synthetic random
     caches (no CAMB run, no model): checks QS_sigma == QS_mu bitwise for arbitrary inputs, checks which bytes of the
     cache each routine writes (struct-layout sanity), and recomputes QS_mu_gen / QS_sigma_gen with an independent
     Python transcription of the general QS formulas.
  A) harness.run_model (harness worker, evol_k) on several models: compare QS_sigma with QS_mu, with QS_sigma_gen, with the
     full-solution Sigma (sigma_eft = -k^2(Phi+Psi)/rho Delta) at sub-horizon k, and, for designer f(R),
     with the analytic QS result Sigma = 1/(1+Omega) built from eft_EFTOmegaV.
Usage: python verify_qs_sigma.py [B|A|all]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REPO = HERE.parents[4]
SCR = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
OUT = SCR / "work/verify/qs_sigma"
TREE = SCR / "trees/A0"

CHILD_B = r"""
import sys, ctypes, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
from camb import eftcamb
from camb.baseconfig import camblib
lib = camblib
names = {"mu": "__eftcamb_abstract_model_MOD_eftcambmodelcomputeqsmu",
         "sigma": "__eftcamb_abstract_model_MOD_eftcambmodelcomputeqssigma",
         "musigma": "__eftcamb_abstract_model_MOD_eftcambmodelcomputeqsmusigma"}
fn = {}
for key, nm in names.items():
    f = getattr(lib, nm)
    f.restype = None
    fn[key] = f
TS = eftcamb.EFTCAMB_timestep_cache
PC = eftcamb.EFTCAMB_parameter_cache
fields = [f for f, _ in TS._fields_]
offs = {f: getattr(TS, f).offset for f in fields}
rng = np.random.default_rng(12345)
fake_self = (ctypes.c_void_p * 4)()  # class container {_data, _vptr}: unused by the routines
out = []
changed_sets = {}
for trial in range(200):
    pc = (PC * 1)()[0]
    c = (TS * 1)()[0]
    pc.h0_Mpc = rng.uniform(1e-4, 4e-4)
    a = rng.uniform(0.1, 1.0)
    k = 10 ** rng.uniform(-3, 0)
    vals = dict(EFTOmegaV=rng.uniform(-0.3, 0.3), EFTOmegaP=rng.uniform(-1, 1), EFTc=rng.uniform(-1e-7, 1e-7),
                adotoa=rng.uniform(1e-4, 1e-3), EFTGamma2V=rng.uniform(-0.3, 0.3), EFTGamma2P=rng.uniform(-1, 1),
                EFTGamma3V=rng.uniform(-0.3, 0.3), EFTGamma4V=rng.uniform(-0.3, 0.3), EFTGamma4P=rng.uniform(-1, 1),
                EFTGamma5V=rng.uniform(-0.3, 0.3), EFTGamma5P=rng.uniform(-1, 1), EFTGamma6V=rng.uniform(-0.3, 0.3))
    if trial % 4 == 0:   # Omega-only (f(R)-like) subset
        for g in ("EFTGamma2V", "EFTGamma2P", "EFTGamma3V", "EFTGamma4V", "EFTGamma4P", "EFTGamma5V", "EFTGamma5P",
                  "EFTGamma6V"):
            vals[g] = 0.0
    if trial % 4 == 1:   # Gamma3=Gamma4=Gamma6=0 (subset where the closed-form QS_mu applies)
        for g in ("EFTGamma3V", "EFTGamma4V", "EFTGamma4P", "EFTGamma6V"):
            vals[g] = 0.0
    vals["Hdot"] = vals["adotoa"] ** 2 * rng.uniform(-1, 0.5)
    vals["Hdotdot"] = vals["adotoa"] ** 3 * rng.uniform(-2, 2)
    for f, v in vals.items():
        setattr(c, f, v)
    before = bytes(c)
    A = ctypes.c_double(a); K = ctypes.c_double(k)
    fn["mu"](ctypes.byref(fake_self), ctypes.byref(A), ctypes.byref(K), ctypes.byref(pc), ctypes.byref(c))
    after_mu = bytes(c)
    fn["sigma"](ctypes.byref(fake_self), ctypes.byref(A), ctypes.byref(K), ctypes.byref(pc), ctypes.byref(c))
    after_sig = bytes(c)
    fn["musigma"](ctypes.byref(fake_self), ctypes.byref(A), ctypes.byref(K), ctypes.byref(pc), ctypes.byref(c))
    after_all = bytes(c)
    def diff(x, y):
        bx = np.frombuffer(x, dtype=np.uint8); by = np.frombuffer(y, dtype=np.uint8)
        idx = np.nonzero(bx != by)[0]
        return sorted({f for f in fields if any(offs[f] <= i < offs[f] + 8 for i in idx)})
    changed_sets.setdefault("mu", set()).update(diff(before, after_mu))
    changed_sets.setdefault("sigma", set()).update(diff(after_mu, after_sig))
    changed_sets.setdefault("musigma", set()).update(diff(after_sig, after_all))
    rec = dict(trial=trial, a=a, k=k, h0=pc.h0_Mpc, vals=vals, QS_mu=c.QS_mu, QS_sigma=c.QS_sigma,
               QS_mu_gen=c.QS_mu_gen, QS_sigma_gen=c.QS_sigma_gen)
    out.append(rec)
print("JSON" + json.dumps(dict(records=out, changed={k: sorted(v) for k, v in changed_sets.items()})))
"""


def gen_formulas(a, k, h0, v):
    """Independent Python transcription of EFTCAMBModelComputeQSMuSigma (06_abstract_EFTCAMB_model.f90:1013-1061)."""
    OV, OP, c, H = v["EFTOmegaV"], v["EFTOmegaP"], v["EFTc"], v["adotoa"]
    Hd, Hdd = v["Hdot"], v["Hdotdot"]
    G2V, G2P, G3V, G4V, G4P = v["EFTGamma2V"], v["EFTGamma2P"], v["EFTGamma3V"], v["EFTGamma4V"], v["EFTGamma4P"]
    G5V, G5P, G6V = v["EFTGamma5V"], v["EFTGamma5P"], v["EFTGamma6V"]
    R0dot = (6 / a**3) * (Hdd - 2 * H * H**2)
    A1 = 2 * (1 + OV) + 4 * G5V
    A2 = -OP * H - h0 * G2V + 2 * H * G4V / a + 4 * H * G5V / a
    A3 = -8 * G6V
    B1 = -1 - 2 * G5V / (1 + OV)
    B2 = 1.0
    B3 = (-OP * H + H * (G4V / a + G4P)) / (1 + OV)
    C1 = OP * H + 2 * H * (G5V / a + G5P)
    C2 = -OP * H / 2 - G2V * h0 / 2 - (1 / (2 * a)) * H * (3 * G3V + G4V - 4 * G5V)
    C3 = (c / a**2 - h0 * H * (G2V / a + G2P) / 2 + G3V * (k**2 / 2 - 3 * (Hd - H**2)) / a**2
          + G4V * (k**2 / 2 - (Hd - H**2)) / a**2 + 2 * (Hd * G5V / a**2 + H**2 * G5P / a))
    M2 = (OP * H * R0dot / 4 - 3 * c * (Hd - H**2) / a**4
          + 1.5 * (h0 / a**2) * ((Hd - H**2) * H * G2P + (Hdd - H * Hd - H**3) * G2V / a)
          + 1.5 * ((Hd - H**2) ** 2 / a**4) * (G4V + 3 * G3V))
    r = a**2 / k**2
    num = B2 * C3 - C1 * B3 + r * B2 * M2
    mu = 2 * num / (A1 * (B3 * C2 - B1 * C3 - B1 * M2 * r) + A2 * (B1 * C1 - B2 * C2) + A3 * (B2 * C3 - B3 * C1 + B2 * M2 * r))
    gam = (B3 * C2 - B1 * C3 - r * B1 * M2) / num
    return mu, 0.5 * mu * (1 + gam)


def part_b() -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        r = subprocess.run([sys.executable, "-c", CHILD_B, str(TREE)], capture_output=True, text=True, env=env,
                           cwd=OUT, timeout=600)
    js = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
    if not js:
        print("child B failed", r.returncode, r.stdout[-2000:], r.stderr[-2000:])
        return {"status": "fail"}
    d = json.loads(js[-1][4:])
    recs = d["records"]
    mu = np.array([x["QS_mu"] for x in recs])
    sg = np.array([x["QS_sigma"] for x in recs])
    mug = np.array([x["QS_mu_gen"] for x in recs])
    sgg = np.array([x["QS_sigma_gen"] for x in recs])
    py = np.array([gen_formulas(x["a"], x["k"], x["h0"], x["vals"]) for x in recs])
    relmug = np.abs(py[:, 0] - mug) / np.maximum(np.abs(mug), 1e-300)
    relsgg = np.abs(py[:, 1] - sgg) / np.maximum(np.abs(sgg), 1e-300)
    n_bitident = int(np.sum(mu == sg))
    gen_diff = np.abs(sgg - mug) / np.maximum(np.abs(mug), 1e-300)
    res = {
        "n_trials": len(recs),
        "n_QS_sigma_bitwise_equal_QS_mu": n_bitident,
        "fields_written_by": d["changed"],
        "py_vs_fortran_mu_gen_maxrel": float(np.max(relmug)),
        "py_vs_fortran_sigma_gen_maxrel": float(np.max(relsgg)),
        "py_vs_fortran_gen_frac_below_1e-10": float(np.mean((relmug < 1e-10) & (relsgg < 1e-10))),
        "gen_sigma_vs_mu_median_reldiff": float(np.median(gen_diff)),
        "gen_sigma_vs_mu_frac_differ_gt_1e-6": float(np.mean(gen_diff > 1e-6)),
    }
    # subsets
    om = [i for i, x in enumerate(recs) if x["trial"] % 4 == 0]
    res["omega_only_subset"] = {
        "n": len(om),
        "QS_sigma_eq_QS_mu": int(np.sum(mu[om] == sg[om])),
        "median_rel_QS_sigma_gen_vs_QS_mu_gen": float(np.median(gen_diff[om])),
    }
    print(json.dumps(res, indent=1))
    return res


MODELS_A = [
    "4_Des_F_R_LCDM_1",
    "4_Des_F_R_wCDM_1",
    "2_PEFT_Omega_const_1",
    "2_PEFT_gamma5_const_1",
    "5_quint_exponential",
]
A_PTS = np.array([0.3, 0.5, 0.8, 1.0])


def part_a() -> dict:
    inis = [REPO / "fortran/eftcamb_test/parameters" / f"{m}.ini" for m in MODELS_A]
    todo = [i for i in inis if not (OUT / f"{i.stem}.npz").exists()]
    if todo:
        st = harness.run_batch(TREE, todo, OUT, jobs=4, evol_k=[0.01, 0.1, 1.0], threads=1, timeout=900)
        for s in st:
            print(Path(s["model_ini"]).stem, s["status"], s["elapsed"])
    res = {}
    for m in MODELS_A:
        f = OUT / f"{m}.npz"
        if not f.exists():
            res[m] = "missing"
            continue
        d = harness.load(f)
        ea = d["evol_a"][0]  # evol_a/evol_k are the cache fields a, k: shape (nk, na)
        idx = [int(np.argmin(np.abs(ea - x))) for x in A_PTS]
        g = lambda name: d[f"evol_{name}"][:, idx]  # noqa: E731
        qm, qs, qmg, qsg, sfull, mfull = (g("QS_mu"), g("QS_sigma"), g("QS_mu_gen"), g("QS_sigma_gen"),
                                          g("sigma_eft"), g("mu"))
        ov = np.interp(ea[idx], d["eft_a_grid"], d["eft_EFTOmegaV"])
        r = {
            "a_used": ea[idx].round(4).tolist(),
            "k": d["evol_k"][:, 0].tolist(),
            "max_abs_QS_sigma_minus_QS_mu": float(np.max(np.abs(qs - qm))),
            "QS_sigma_k1": qs[2].round(5).tolist(),
            "QS_sigma_gen_k1": qsg[2].round(5).tolist(),
            "full_sigma_eft_k1": sfull[2].round(5).tolist(),
            "full_mu_k1": mfull[2].round(5).tolist(),
            "QS_mu_k1": qm[2].round(5).tolist(),
            "inv_1_plus_OmegaV": (1 / (1 + ov)).round(5).tolist(),
            "maxrel_QS_sigma_vs_full_k1": float(np.max(np.abs(qs[2] / sfull[2] - 1))),
            "maxrel_QS_sigma_gen_vs_full_k1": float(np.max(np.abs(qsg[2] / sfull[2] - 1))),
            "maxrel_QS_mu_vs_full_mu_k1": float(np.max(np.abs(qm[2] / mfull[2] - 1))),
        }
        res[m] = r
        print(m, json.dumps(r))
    return res


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    allres = {}
    if which in ("B", "all"):
        allres["B"] = part_b()
    if which in ("A", "all"):
        allres["A"] = part_a()
    (OUT / f"verify_qs_sigma_{which}.json").write_text(json.dumps(allres, indent=1))

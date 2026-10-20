"""Adversarial verification (steplog-d4-and-integral): model-level impact, A0 vs a copy of A0 with ONLY the steplog
d4 + integral fix (work/verify/steplog_phys/fixtree, built by steplog_build_fixtree.py).

Runs (threads=1, <= 4 concurrent, all via harness.run_model):
  C0  control: 4_Des_F_R_wCDM_1 unmodified (no steplog) -> must be bit-identical between trees
  P1  Pure EFT Omega = steplog(0, 1.5, 0.5, 0.5), mass stability F         (finding's case; NaN onset a = 0.642)
  P2  Pure EFT Omega = steplog(0, 0.1, 0.6, 0.6), mass stability F         (moderate; NaN onset a = 0.860)
  P3  P2 with EFTCAMB_skip_stability = T: does the NaN reach the perturbations / Cls?
  P4  Pure EFT Omega = steplog(0, 0.1, 0.8, 0.6), mass stability F         (no NaN: aT e^delta^2 = 1.15 > 1)
  P5  P4 with EFT_mass_stability = T (mass stability uses EFTcdotdot, which contains EFTOmegaPPPP)
  R1  RPH (3_RPH_MB_const_1) Planck mass = steplog(0, 0.1, 0.6, 0.6), mass stability F
  R2  RPH (3_RPH_MK_const_1) Planck mass = steplog(0, 0.1, 0.6, 0.6); R3 same with (0, 0.5, 0.3, 0.6)
  R4  RPH (3_RPH_KT_power_1) tensor alpha_T = steplog(0, -0.01, 0.6, 0.6)
  D1  designer f(R) (4_Des_F_R_wCDM_1), w = steplog(-1, -0.8, 0.3, 0.5)  (finding's case)
  D2  designer f(R), w = steplog(-1.0, -0.9, 0.6, 0.3), EFTB0 = 0.01
  S1  ShiftSym alphaB (audit ini 007p7-shiftsym.ini), w = steplog(-1.0, -0.9, 0.5, 0.4)
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "verify" / "steplog_phys"
FIX = WORK / "fixtree"
P = harness.PARAMS_DIR
INIS = harness.REPO / "fortran" / "eftcamb_test" / "audit" / "inis"


def om(v1, v2, at, de, pre="EFTOmega"):
    return {f"{pre}_v1": str(v1), f"{pre}_v2": str(v2), f"{pre}_at": str(at), f"{pre}_delta": str(de)}


RUNS = {
    "C0": (P / "4_Des_F_R_wCDM_1.ini", {}),
    "P1": (P / "2_PEFT_Omega_lin_1.ini", {"EFT_mass_stability": "F", "PureEFTmodelOmega": "8", **om(0.0, 1.5, 0.5, 0.5)}),
    "P2": (P / "2_PEFT_Omega_lin_1.ini", {"EFT_mass_stability": "F", "PureEFTmodelOmega": "8", **om(0.0, 0.1, 0.6, 0.6)}),
    "P3": (P / "2_PEFT_Omega_lin_1.ini", {"EFT_mass_stability": "F", "PureEFTmodelOmega": "8", **om(0.0, 0.1, 0.6, 0.6),
                                          "EFTCAMB_skip_stability": "T"}),
    "P4": (P / "2_PEFT_Omega_lin_1.ini", {"EFT_mass_stability": "F", "PureEFTmodelOmega": "8", **om(0.0, 0.1, 0.8, 0.6)}),
    "P5": (P / "2_PEFT_Omega_lin_1.ini", {"EFT_mass_stability": "T", "PureEFTmodelOmega": "8", **om(0.0, 0.1, 0.8, 0.6)}),
    "R1": (P / "3_RPH_MB_const_1.ini", {"EFT_mass_stability": "F", "RPHmassPmodel": "8",
                                        **om(0.0, 0.1, 0.6, 0.6, pre="RPHmassP")}),
    "R2": (P / "3_RPH_MK_const_1.ini", {"EFT_mass_stability": "F", "RPHmassPmodel": "8",
                                        **om(0.0, 0.1, 0.6, 0.6, pre="RPHmassP")}),
    "R3": (P / "3_RPH_MK_const_1.ini", {"EFT_mass_stability": "F", "RPHmassPmodel": "8",
                                        **om(0.0, 0.5, 0.3, 0.6, pre="RPHmassP")}),
    "R4": (P / "3_RPH_KT_power_1.ini", {"EFT_mass_stability": "F", "RPHtensormodel": "8",
                                        **om(0.0, -0.01, 0.6, 0.6, pre="RPHtensor")}),
    "D1": (P / "4_Des_F_R_wCDM_1.ini", {"EFTwDE": "8", **om(-1.0, -0.8, 0.3, 0.5, pre="EFTw")}),
    "D2": (P / "4_Des_F_R_wCDM_1.ini", {"EFTwDE": "8", "EFTB0": "0.01", **om(-1.0, -0.9, 0.6, 0.3, pre="EFTw")}),
    "S1": (INIS / "007p7-shiftsym.ini", {"EFTwDE": "8", **om(-1.0, -0.9, 0.5, 0.4, pre="EFTw")}),
}


def job(args):
    name, tag, tree = args
    ini, ov = RUNS[name]
    npz = WORK / f"{name}_{tag}.npz"
    st = harness.run_model(tree, ini, npz, overrides=ov, threads=1, timeout=1500)
    log = Path(st["log"]).read_text(errors="replace").splitlines() if st.get("log") else []
    msg = [t.strip() for t in log if ("Instability detected" in t or "Model has Nan" in t or "theory unstable" in t
                                      or "unstable" in t.lower())][:3]
    first_nan = next((t.strip() for t in log if "is NaN" in t), "")
    return name, tag, {"status": st["status"], "msg": msg, "first_nan": first_nan, "npz": str(npz)}


def main():
    only = sys.argv[1].split(",") if len(sys.argv) > 1 else list(RUNS)
    jobs = [(n, t, tr) for n in only for t, tr in (("A0", "A0"), ("FIX", FIX))]
    res = {}
    with ThreadPoolExecutor(4) as ex:
        for name, tag, r in ex.map(job, jobs):
            res.setdefault(name, {})[tag] = r
            print(name, tag, r["status"], r["first_nan"], r["msg"], flush=True)
    summary = {}
    for name in only:
        a, f = res[name]["A0"], res[name]["FIX"]
        s = {"A0": {k: a[k] for k in ("status", "first_nan", "msg")}, "FIX": {k: f[k] for k in ("status", "first_nan", "msg")}}
        if a["status"] == "ok" and f["status"] == "ok":
            da, df = harness.load(a["npz"]), harness.load(f["npz"])
            c = harness.compare(da, df)
            s["bit_identical"] = c["bit_identical"]
            s["worst_cl"] = c["worst_cl"]
            s["worst_pk"] = c["worst_pk"]
            s["sigma8"] = [float(np.atleast_1d(da["sigma8"])[-1]), float(np.atleast_1d(df["sigma8"])[-1])]
            s["bg_H_rel"] = c.get("bg_H")
            s["bg_tau_rel"] = c.get("bg_tau")
            s["cl_spectra"] = c["cl"]
            ag = da["eft_a_grid"]
            for fld in ("eft_grhov_t", "eft_adotoa", "eft_EFTOmegaPPPP", "eft_EFTcdotdot", "eft_EFTOmegaV"):
                if fld in da and fld in df:
                    s[fld + "@a=1"] = [float(da[fld][-1]), float(df[fld][-1])]
                    s[fld + "@a=" + f"{ag[-2]:.4f}"] = [float(da[fld][-2]), float(df[fld][-2])]
            eftd = {k: v for k, v in c.get("eft", {}).items() if v and v > 0}
            s["eft_maxreldiff_top"] = dict(sorted(eftd.items(), key=lambda kv: -kv[1])[:6])
            meta_a, meta_f = da["meta"], df["meta"]
            s["derived_A0"] = {k: meta_a.get("derived", {}).get(k) for k in ("age", "rstar", "thetastar", "zdrag")}
            s["derived_FIX"] = {k: meta_f.get("derived", {}).get(k) for k in ("age", "rstar", "thetastar", "zdrag")}
        summary[name] = s
    out = WORK / ("models_" + "_".join(only) + ".json")
    out.write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()

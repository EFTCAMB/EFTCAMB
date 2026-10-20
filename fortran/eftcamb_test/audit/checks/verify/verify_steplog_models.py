"""Model-level verification of steplog-d4-and-integral through the harness (CAMB runs, threads=1).

 C1  Pure EFT, Omega = steplog(0, 0.05, 0.6, 0.3): d4 NaN for a > 0.6 e^0.09 = 0.656 -> expect rejection near a = 0.656
 C2  control: same with delta = 0.75 (aT e^delta^2 = 1.053 > 1): expect ok; compare eft_EFTOmegaPPPP with the exact d4
     and with a finite difference of the run's own eft_EFTOmegaPPP
 C3  control C2 with EFT_mass_stability = T (mass stability uses EFTcdotdot, which contains EFTOmegaPPPP)
 C1b Pure EFT, Omega = steplog(0, 0.05, 0.45, 0.75): d4 NaN for a > 0.79 (otherwise like control C2)
 C5  the t4a-param repro case steplog(0, 1.5, 0.5, 0.5)
 C4  designer f(R) with EFTwDE = 8, w = steplog(-1, -0.9, 0.5, 0.3): grhov_t(a=1)/grhov_t(a<1) vs exact ratio of I(a);
     CAMB's own H(z=0) vs H0 input

usage: python verify_steplog_models.py [--tree A0|path] [--tag A0]
"""

import argparse
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
           "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/steplog")


def wfun(v1, v2, at, de):
    v1, v2, at, de = (mp.mpf(str(t)) for t in (v1, v2, at, de))

    def w(a):
        u = (mp.log(a) - mp.log(at)) / de
        return (v2 - v1) / 2 * u / mp.sqrt(1 + u * u) + (v1 + v2) / 2
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default="A0")
    ap.add_argument("--tag", default="A0")
    ap.add_argument("--only", default="C1,C2,C3,C4")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    mp.mp.dps = 30
    peft = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"
    dfr = harness.PARAMS_DIR / "4_Des_F_R_wCDM_1.ini"
    runs = {
        "C1": (peft, {"EFT_mass_stability": "F", "PureEFTmodelOmega": "8", "EFTOmega_v1": "0.0",
                      "EFTOmega_v2": "0.05", "EFTOmega_at": "0.6", "EFTOmega_delta": "0.3"}),
        "C2": (peft, {"EFT_mass_stability": "F", "PureEFTmodelOmega": "8", "EFTOmega_v1": "0.0",
                      "EFTOmega_v2": "0.05", "EFTOmega_at": "0.6", "EFTOmega_delta": "0.75"}),
        "C3": (peft, {"EFT_mass_stability": "T", "PureEFTmodelOmega": "8", "EFTOmega_v1": "0.0",
                      "EFTOmega_v2": "0.05", "EFTOmega_at": "0.6", "EFTOmega_delta": "0.75"}),
        "C1b": (peft, {"EFT_mass_stability": "F", "PureEFTmodelOmega": "8", "EFTOmega_v1": "0.0",
                       "EFTOmega_v2": "0.05", "EFTOmega_at": "0.45", "EFTOmega_delta": "0.75"}),
        "C5": (peft, {"EFT_mass_stability": "F", "PureEFTmodelOmega": "8", "EFTOmega_v1": "0.0",
                      "EFTOmega_v2": "1.5", "EFTOmega_at": "0.5", "EFTOmega_delta": "0.5"}),
        "C4": (dfr, {"EFTwDE": "8", "EFTw_v1": "-1.0", "EFTw_v2": "-0.9", "EFTw_at": "0.5", "EFTw_delta": "0.3"}),
    }
    summary = {}
    for name in args.only.split(","):
        ini, ov = runs[name]
        npz = OUT / f"{name}_{args.tag}.npz"
        st = harness.run_model(args.tree, ini, npz, overrides=ov, threads=1)
        log = Path(st["log"]).read_text(errors="replace").splitlines() if st.get("log") else []
        inst = [t.strip() for t in log if "Instability detected" in t or "Model has Nan" in t
                or "is NaN" in t][:2]
        s = {"status": st["status"], "log": inst}
        if st["status"] == "ok":
            d = harness.load(npz)
            a = d["eft_a_grid"]
            if name in ("C2", "C3"):
                w = wfun(0.0, 0.05, 0.6, 0.75)
                sel = [i for i in range(len(a)) if a[i] > 0.02]
                pppp, ppp = d["eft_EFTOmegaPPPP"], d["eft_EFTOmegaPPP"]
                fd = np.gradient(ppp, a)
                rows = []
                for i in sel[:: max(1, len(sel) // 8)] + [len(a) - 1]:
                    ex = float(mp.diff(w, mp.mpf(a[i]), 4))
                    rows.append((float(a[i]), float(pppp[i]), ex, float(fd[i])))
                s["d4_rows(a,obs,exact,fd_of_run_PPP)"] = rows
                s["max_rel_err_obs_vs_exact"] = max(abs(r[1] - r[2]) / abs(r[2]) for r in rows)
                s["cl_TT_2000"] = float(d["cl_total"][2000, 0])
            if name == "C4":
                w = wfun(-1.0, -0.9, 0.5, 0.3)

                def integ(x):
                    x = mp.mpf(x)
                    return x**2 * mp.exp(3 * mp.quad(lambda n: 1 + w(mp.exp(n)), [mp.log(x), 0]))
                g = d["eft_grhov_t"]
                s["grhov_t_last3"] = [float(t) for t in g[-3:]]
                s["a_last3"] = [float(t) for t in a[-3:]]
                s["ratio_obs_g(1)/g(a_-1)"] = float(g[-1] / g[-2])
                s["ratio_exact_I(1)/I(a_-1)"] = float(1 / integ(a[-2]))
                s["ratio_obs_g(a_-1)/g(a_-2)"] = float(g[-2] / g[-3])
                s["ratio_exact_I(a_-1)/I(a_-2)"] = float(integ(a[-2]) / integ(a[-3]))
                s["bg_z0"] = float(d["bg_z"][0])
                s["bg_H_z0"] = float(d["bg_H"][0])
                s["H0_input"] = 70.0
                s["adotoa_last2"] = [float(t) for t in d["eft_adotoa"][-2:]]
                s["sigma8"] = float(np.atleast_1d(d["sigma8"])[-1])
        summary[name] = s
        print(name, json.dumps(s, indent=1))
    (OUT / f"models_{args.tag}.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

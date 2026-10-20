"""t8b check: effect of the cmbmain.f90 GetTransfer EFTCAMB MOD (transfer_high_precision no longer tightens atol by
1/100 when EFTflag /= 0). GetTransfer is only used for transfer-only k values, so run with WantCls=False,
WantTransfer=True, kmax=10/Mpc, high_precision=T, and compare P(k) of A0 (code under audit) and A0_hp (upstream
atol/100 restored) with a boosted reference (A0_hp, AccuracyBoost=2).
Usage: python hp_transfer_only.py  (runs the three configurations and prints max |dP/P|)"""

import sys

import numpy as np


def main(args):
    import camb

    from t8b_common import read_model

    out = {}
    for ini in args["inis"]:
        try:
            cp = read_model(camb, ini)
            cp.WantCls = False
            cp.SourceWindows = []  # upstream SIGSEGV with windows and WantCls=False (t10a)
            cp.DoLensing = False
            cp.Transfer.high_precision = True
            cp.set_matter_power(redshifts=[1.0, 0.0], kmax=10.0, k_per_logint=0)
            cp.NonLinear = camb.model.NonLinear_none
            if args.get("boost"):
                cp.Accuracy.AccuracyBoost = args["boost"]
            res = camb.get_results(cp)
            kh, z, pk = res.get_matter_power_spectrum(minkh=1e-3, maxkh=6.0, npoints=300)
            out[ini] = {"kh": kh.tolist(), "pk": pk.tolist(), "sigma8": res.get_sigma8().tolist()}
        except Exception as e:  # noqa: BLE001
            out[ini] = {"error": repr(e)[:200]}
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, WORK, run_in_tree

    inis = sys.argv[1:] or ["1_EFT_GR.ini", "2_PEFT_gamma1_const_1.ini", "4_Des_F_R_wCDM_1.ini", "5_Kmouflage_1.ini"]
    hp = str(WORK / "A0_hp")
    runs = {"A0": ("A0", None), "A0_hp": (hp, None), "ref": (hp, 2.0)}
    res = {}
    for name, (tree, boost) in runs.items():
        res[name] = run_in_tree(tree, HERE / "hp_transfer_only.py", {"inis": inis, "boost": boost}, timeout=3000)
    (WORK / "out" / "hp_transfer_only.json").write_text(json.dumps(res))
    for ini in inis:
        line = f"{ini:28s}"
        try:
            ref = np.array(res["ref"][ini]["pk"])
            for name in ("A0", "A0_hp"):
                p = np.array(res[name][ini]["pk"])
                line += f" {name}: max|dP/P|={np.max(np.abs(p / ref - 1)):.2e}"
            dA = np.max(np.abs(np.array(res["A0"][ini]["pk"]) / np.array(res["A0_hp"][ini]["pk"]) - 1))
            line += f"  A0 vs A0_hp={dA:.2e}"
            s = [res[n][ini]["sigma8"][-1] for n in ("A0", "A0_hp", "ref")]
            line += f"  sigma8 A0/ref-1={s[0] / s[2] - 1:+.1e} hp/ref-1={s[1] / s[2] - 1:+.1e}"
        except Exception as e:  # noqa: BLE001
            line += " " + repr(e)[:80] + " " + str({n: res[n].get(ini, res[n]).get("error") for n in res})[:300]
        print(line)

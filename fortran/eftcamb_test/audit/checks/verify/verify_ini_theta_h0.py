"""Independent check of ini-theta-h0-solved-without-eft.

1) Write ini with thetastar via DEFAULT(model.ini) (no flattening), read via camb.read_ini -> H0_ini.
2) Same ini with hubble=H0_ini and EFTCAMB_flag=0 -> theta must equal request (H0 solved for GR background).
3) Same ini with hubble=H0_ini and EFT on -> theta off (as reported).
4) Brute-force: scipy brentq on H0 using hubble= ini runs with EFT on -> correct H0.
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
from scipy.optimize import brentq
TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")
OUT = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
       "scratchpad/audit/work/verify/ini-theta")
PAR = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"
sys.path.insert(0, TREE)
with harness.slot(1):
    import camb
    TH = 0.010411
    def theta_from(model, extra):
        lines = [f"DEFAULT({PAR}/{model}.ini)"] + [f"{k} = {v}" for k, v in extra.items()]
        base = {"get_scalar_cls": "F", "get_tensor_cls": "F", "get_transfer": "F", "num_redshiftwindows": "0",
                "feedback_level": "0", "EFTCAMB_write_background": "F"}
        lines += [f"{k} = {v}" for k, v in base.items()]
        p = f"{OUT}/{model}_{abs(hash(tuple(extra.items())))}.ini"
        open(p, "w").write("\n".join(lines) + "\n")
        # hubble key must be removed when thetastar set: CAMB errors if both; use hubble= empty not possible,
        # so the thetastar ini is written without DEFAULT hubble via a flattened copy below.
        cp = camb.read_ini(p)
        return cp, camb.get_background(cp).get_derived_params()["thetastar"] / 100
    def flat_theta_ini(model):
        k = harness.parse_ini(f"{PAR}/{model}.ini")
        k.pop("hubble", None)
        return k
    for model in ("2_PEFT_wDE_wCDM_1", "5_quint_exponential", "5_Kmouflage_1", "4_Des_F_R_wCDM_1"):
        k = flat_theta_ini(model)
        k |= {"thetastar": str(TH), "get_scalar_cls": "F", "get_tensor_cls": "F", "get_transfer": "F",
              "num_redshiftwindows": "0", "feedback_level": "0", "EFTCAMB_write_background": "F"}
        p = f"{OUT}/{model}_theta.ini"
        open(p, "w").write("\n".join(f"{a} = {b}" for a, b in k.items()) + "\n")
        cp = camb.read_ini(p)
        H0i = cp.H0
        th_ini = camb.get_background(cp).get_derived_params()["thetastar"] / 100
        _, th_gr = theta_from(model, {"hubble": repr(H0i), "EFTflag": "0"})
        _, th_eft = theta_from(model, {"hubble": repr(H0i)})
        f = lambda h: theta_from(model, {"hubble": repr(h)})[1] - TH
        try:
            Hc = brentq(f, 50, 95, xtol=1e-4)
        except Exception as e:  # noqa: BLE001
            print("brentq failed", repr(e)[:200])
            Hc = float("nan")
        print(f"{model:20s} EFTflag={cp.EFTCAMB.EFTflag} H0_ini={H0i:.4f} theta_ini={th_ini:.7f} "
              f"rel={th_ini/TH-1:+.2e} | hubble=H0_ini,EFTflag=0: rel={th_gr/TH-1:+.2e} | EFT on: rel={th_eft/TH-1:+.2e}"
              f" | brentq H0={Hc:.4f}")

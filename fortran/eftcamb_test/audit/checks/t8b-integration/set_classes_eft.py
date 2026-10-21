"""t8b check: camb/model.py set_classes EFTCAMB MOD. The comment says that when dark_energy_model != 'EFTCAMB'
EFTCAMB is initialised with EFTflag = 0, but the elif branch leaves an active EFTCAMB untouched, so a later
set_dark_energy / set_classes(dark_energy_model='ppf') is silently ignored. Usage: python set_classes_eft.py TREE"""

import sys


def main(args):
    import camb

    from t8b_common import read_model

    cp = read_model(camb, "2_PEFT_wDE_wCDM_1.ini")
    out = {"EFTflag_before": int(cp.EFTCAMB.EFTflag)}
    H_eft = float(camb.get_background(cp, no_thermo=True).hubble_parameter(1.0))
    cp.set_dark_energy(w=-0.7, dark_energy_model="ppf")
    out["EFTflag_after_set_dark_energy_ppf"] = int(cp.EFTCAMB.EFTflag)
    out["DarkEnergy_class"] = type(cp.DarkEnergy).__name__
    H_after = float(camb.get_background(cp, no_thermo=True).hubble_parameter(1.0))
    cp2 = camb.set_params(H0=cp.H0, ombh2=cp.ombh2, omch2=cp.omch2, omk=cp.omk, mnu=0.0, num_massive_neutrinos=0,
                          nnu=cp.N_eff, w=-0.7, dark_energy_model="ppf")
    H_w07 = float(camb.get_background(cp2, no_thermo=True).hubble_parameter(1.0))
    out.update({"H(z=1) EFT model": H_eft, "H(z=1) after set_dark_energy(w=-0.7)": H_after, "H(z=1) plain w=-0.7": H_w07})
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, run_in_tree

    print(json.dumps(run_in_tree(sys.argv[1], HERE / "set_classes_eft.py", {}), indent=1))

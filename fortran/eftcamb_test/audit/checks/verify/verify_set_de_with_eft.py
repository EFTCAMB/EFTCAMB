"""Independent check: EFT initialised from dict, then set_dark_energy_w_a / set_classes('EFTCAMB') ; compare TT."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import harness
os.environ["OMP_NUM_THREADS"] = "1"
TREE = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0"
sys.path.insert(0, TREE)
with harness.slot(1):
    import camb, numpy as np
    def run(cp):
        cp.set_for_lmax(1500)
        r = camb.get_results(cp)
        return r.hubble_parameter(1.0), r.get_cmb_power_spectra(cp, CMB_unit="muK")["total"][2:1500, 0]
    base = dict(H0=67.0, ombh2=0.022, omch2=0.12, As=2e-9, ns=0.965)
    eft = {"EFTflag": 1, "PureEFTmodel": 1, "EFTwDE": 1, "EFTw0": -0.9, "PureEFTmodelOmega": 0,
           "PureEFTmodelGamma1": 0, "PureEFTmodelGamma2": 0, "PureEFTmodelGamma3": 0,
           "PureEFTmodelGamma4": 0, "PureEFTmodelGamma5": 0, "PureEFTmodelGamma6": 0, "PureEFTHorndeski": False,
           "EFTCAMB_feedback_level": 0}
    cp = camb.set_params(**base)
    cp.EFTCAMB.initialize_parameters(cp, eft, print_header=False)
    h_eft, tt_eft = run(cp)
    a = np.linspace(1e-4, 1, 50); w = -0.7 + 0 * a
    cp.set_dark_energy_w_a(a, w, dark_energy_model="ppf")
    h_sw, tt_sw = run(cp)
    flag_after_wa = cp.EFTCAMB.EFTflag
    ref = camb.set_params(**base, w=-0.7, dark_energy_model="ppf")
    h_ref, tt_ref = run(ref)
    cpw = camb.set_params(**base, w=-0.9, dark_energy_model="fluid")
    h_w9, tt_w9 = run(cpw)
    cp2 = camb.set_params(**base)
    cp2.EFTCAMB.initialize_parameters(cp2, eft, print_header=False)
    cp2.set_classes(dark_energy_model="EFTCAMB")
    flag_after_eftcamb = cp2.EFTCAMB.EFTflag
    f = lambda x, y: np.max(np.abs(x / y - 1))
    print(f"EFT wCDM(-0.9): EFTflag=1 H(z=1)={h_eft:.4f}; plain fluid w=-0.9 H={h_w9:.4f}, max|dTT/TT|={f(tt_eft, tt_w9):.2e}")
    print(f"after set_dark_energy_w_a(ppf,-0.7): EFTflag={flag_after_wa} H={h_sw:.4f}, max|dTT| vs EFT={f(tt_sw, tt_eft):.2e}")
    print(f"plain ppf w=-0.7: H={h_ref:.4f}; switched vs ref max|dTT/TT|={f(tt_sw, tt_ref):.2e}")
    print(f"set_classes('EFTCAMB') on active EFT -> EFTflag={flag_after_eftcamb}")

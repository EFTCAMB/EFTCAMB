"""Verifier: 21cm + EFT through (a) run_ini (driver), (b) pure-Python SourceWindows path, (c) EFT GR-limit vs GR."""
import os, sys, subprocess
from pathlib import Path
os.environ["OMP_NUM_THREADS"] = "1"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness
TREE = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0"
OUT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/21cm")
PARAMS = HERE.parents[2] / "parameters"
mode = sys.argv[1] if len(sys.argv) > 1 else "main"
sys.path.insert(0, TREE)
if mode == "runini":
    import camb
    ini = OUT / "runini.ini"
    ini.write_text(f"DEFAULT({PARAMS/'2_PEFT_Omega_const_1.ini'})\nget_tensor_cls = F\nget_transfer = F\nfeedback_level = 0\n"
                   f"EFTCAMB_write_background = F\nDo21cm = T\nnum_redshiftwindows = 1\nredshift(1) = 1.0\n"
                   f"redshift_kind(1) = 21cm\nredshift_sigma_Mhz(1) = 5\nl_max_scalar = 600\noutput_root = {OUT}/runini\n")
    camb.run_ini(str(ini), no_validate=True)
    print("run_ini returned normally")
    sys.exit(0)
with harness.slot(1):
    if mode == "main":
        r = subprocess.run([sys.executable, __file__, "runini"], capture_output=True, text=True, cwd=OUT, timeout=600)
        print("run_ini subprocess rc=", r.returncode, "| tail:", (r.stdout + r.stderr).strip().splitlines()[-3:])
    import camb
    from camb.sources import GaussianSourceWindow
    def cl21(eft, om0=None):
        p = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965, lmax=600)
        p.SourceWindows = [GaussianSourceWindow(redshift=1.0, source_type="21cm", sigma=0.02)]
        p.SourceTerms.line_distortions = True
        p.SourceTerms.limber_windows = True
        if eft:
            p.EFTCAMB.initialize_parameters(p, dict(EFTflag=1, PureEFTmodel=1, PureEFTmodelOmega=1, EFTOmega0=om0,
                EFTOmegaExp=1.0, EFTwDE=0, EFTCAMB_feedback_level=0, PureEFTHorndeski=False), print_header=False)
        try:
            res = camb.get_results(p)
        except Exception as e:  # noqa
            return f"raised {type(e).__name__}: {e}"
        return res.get_source_cls_dict(raw_cl=True)["W1xW1"]
    gr = cl21(False)
    for om0 in (1e-8, 0.1, 1.0):
        c = cl21(True, om0)
        if isinstance(c, str):
            print("Omega0", om0, c)
        else:
            print(f"pure-python EFT Omega0={om0}: no error; C_l l=10 {c[10]:.4e} vs GR {gr[10]:.4e}, ratio l=10 {c[10]/gr[10]:.5f} l=100 {c[100]/gr[100]:.5f}")

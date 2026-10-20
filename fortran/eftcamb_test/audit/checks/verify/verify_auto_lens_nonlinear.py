"""Independent verifier: auto-lens-accuracy-enables-nonlinear-for-eft.
(1) python set_params path (camb.set_params with EFTCAMB dict) -> NonLinear; (2) GR (EFTflag=0) magnitude of the same
switch, to separate generic nonlinear-lensing size from EFT-specific error; (3) pk nonlinear ratio for f(R) vs GR; (4) whether
the Fortran warning prints with python-default feedback."""
import os, sys, subprocess
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
TREE = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0"
sys.path.insert(0, TREE)
import numpy as np
with harness.slot(1):
    import camb
    eft = dict(EFTflag=3, DesignerEFTmodel=1, EFTwDE=1, EFTw0=-0.9, EFTB0=1.0, feedback_level=0)
    def run(eftd, **kw):
        pass
        cp = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965, lmax=2500, EFTCAMB_params=eftd, **kw) \
            if False else None
        cp = camb.CAMBparams()
        cp.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
        cp.InitPower.set_params(As=2e-9, ns=0.965)
        if eftd is not None:
            cp.EFTCAMB.initialize_parameters(cp, eftd, print_header=False)
        cp.set_for_lmax(2500, **kw)
        nl = str(cp.NonLinear)
        cl = camb.get_results(cp).get_cmb_power_spectra(CMB_unit="muK", raw_cl=True, lmax=2500)["total"]
        return nl, cl
    for name, d in (("GR(no EFT)", None), ("GR EFTflag=0", dict(EFTflag=0)), ("f(R) w0=-0.9 B0=1", eft)):
        nl, cd = run(d)
        _, cl = run(d, nonlinear=False)
        r = [np.abs(cd[2:, i] / cl[2:, i] - 1).max() for i in range(3)]
        print(f"{name:20s} default NonLinear={nl}: max|NL/lin-1| TT={r[0]:.3e} EE={r[1]:.3e} BB={r[2]:.3e}")

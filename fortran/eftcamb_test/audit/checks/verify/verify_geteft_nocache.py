"""Independent verifier: get_eft_functions / get_scale_evolution on results without an initialised EFT cache."""
import os, subprocess, sys
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness

TREE = sys.argv[1] if len(sys.argv) > 1 else ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
                                              "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
INI = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/2_PEFT_Omega_const_1.ini"
PRE = r"""
import sys, faulthandler; faulthandler.enable()
sys.path.insert(0, sys.argv[1])
import camb, numpy as np
from camb import model
"""
CASES = {
 "ctrl_eft_background": "pars=camb.read_ini(sys.argv[2], no_validate=True)\nres=camb.get_background(pars)\n"
     "f,a=pars.EFTCAMB.get_eft_functions(res,np.array([0.5,1.0]))\nprint('RES adotoa',a['adotoa'],'EFTOmegaV',a['EFTOmegaV'])",
 "ctrl_eft_full": "pars=camb.read_ini(sys.argv[2], no_validate=True)\npars.WantCls=True\nres=camb.get_results(pars)\n"
     "f,a=pars.EFTCAMB.get_eft_functions(res,np.array([0.5,1.0]))\nprint('RES adotoa',a['adotoa'])",
 "gr_full_results": "pars=camb.set_params(H0=67.5,ombh2=0.022,omch2=0.122,As=2e-9,ns=0.965)\nres=camb.get_results(pars)\n"
     "print('RES EFTflag',pars.EFTCAMB.EFTflag)\nf,a=pars.EFTCAMB.get_eft_functions(res,np.array([0.5,1.0]))\nprint('RES adotoa',a['adotoa'])",
 "gr_background_mnu0": "pars=camb.set_params(H0=67.5,ombh2=0.022,omch2=0.122,mnu=0,As=2e-9,ns=0.965)\nres=camb.get_background(pars)\n"
     "f,a=pars.EFTCAMB.get_eft_functions(res,np.array([0.5,1.0]))\nprint('RES adotoa',a['adotoa'])",
 "eft_noout_results": "pars=camb.read_ini(sys.argv[2], no_validate=True)\npars.WantCls=False\npars.WantTransfer=False\npars.DoLensing=False\n"
     "res=camb.get_results(pars)\nprint('RES EFTflag',pars.EFTCAMB.EFTflag)\nf,a=pars.EFTCAMB.get_eft_functions(res,np.array([0.5,1.0]))\nprint('RES adotoa',a['adotoa'])",
 "eft_noout_scale_evol": "pars=camb.read_ini(sys.argv[2], no_validate=True)\npars.WantCls=False\npars.WantTransfer=False\npars.DoLensing=False\n"
     "res=camb.get_results(pars)\nf,a=pars.EFTCAMB.get_scale_evolution(res,np.array([0.01]),np.array([0.5,1.0]))\n"
     "print('RES adotoa',np.ravel(a['adotoa']),'a',np.ravel(a['a']))",
}
env = dict(os.environ, OMP_NUM_THREADS="1"); env.pop("PYTHONPATH", None)
only = sys.argv[2:] or list(CASES)
for name in only:
    with harness.slot(1):
        r = subprocess.run([sys.executable, "-c", PRE + CASES[name], TREE, INI], env=env, capture_output=True,
                           text=True, timeout=500, cwd="/private/tmp")
    res = [l for l in r.stdout.splitlines() if l.startswith("RES")]
    err = [l for l in r.stderr.splitlines() if "Error" in l or "Fatal" in l or "File" in l][:4]
    print(f"{name}: rc={r.returncode} {res} {err}", flush=True)

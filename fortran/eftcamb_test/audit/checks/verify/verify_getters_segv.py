"""Adversarial verification of F-020 (eft-getters-gr-segfault).

Each case runs in its own subprocess (one CAMB slot, OMP_NUM_THREADS=1) importing camb from the chosen tree.
Prints, per tree and case: return code / signal, the last stdout lines and the stderr tail.

Cases:
  default_<getter>      : camb.CAMBparams() (EFTflag=0) -> getter
  setparams_gr_<getter> : camb.set_params(H0=67.5, EFTflag=0) (documented GR path through EFT_do_set) -> getter
  valid_<getter>        : camb.read_ini(2_PEFT_Omega_const_1.ini) (model allocated, control) -> getter
  failmodel_*           : initialize_parameters({'EFTflag':1,'PureEFTmodel':9}) raises -> getter / get_background /
                          get_results
  failomega_*           : initialize_parameters({'EFTflag':1,'PureEFTmodel':1,'PureEFTmodelOmega':99}) raises
                          (model allocated, sub-function missing) -> getter / get_results
  reinit_ok             : failed init, then a valid initialize_parameters on the same object -> get_results works?

Usage: python verify_getters_segv.py OUTFILE TREE [TREE ...] [--cases c1,c2]
"""

import json
import os
import signal
import subprocess
import sys

sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

REPO = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB"
INI = REPO + "/fortran/eftcamb_test/parameters/2_PEFT_Omega_const_1.ini"

CHILD = r"""
import sys, os
tree, case = sys.argv[1], sys.argv[2]
sys.path.insert(0, tree)
os.chdir(os.environ['VERIFY_CWD'])
import camb
def show(p, tag=''):
    print('STATE', tag, 'EFTflag =', p.EFTCAMB.EFTflag, flush=True)
kind, _, what = case.partition('_')
if kind == 'default':
    p = camb.CAMBparams()
elif kind == 'setparamsgr':
    try:
        p = camb.set_params(H0=67.5, EFTflag=0, lmax=300)
    except TypeError:
        print('NOTE set_params TypeError, retry with YHe=0.2454', flush=True)
        p = camb.set_params(H0=67.5, EFTflag=0, lmax=300, YHe=0.2454)
elif kind == 'valid':
    p = camb.read_ini(sys.argv[3])
elif kind in ('failmodel', 'failomega', 'reinit'):
    p = camb.CAMBparams()
    d = {'EFTflag': 1, 'PureEFTmodel': 9} if kind != 'failomega' else \
        {'EFTflag': 1, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 99}
    try:
        p.EFTCAMB.initialize_parameters(p, d, print_header=False)
        print('NO EXCEPTION', flush=True)
    except Exception as e:
        print('initialize_parameters raised', type(e).__name__, flush=True)
    if kind == 'reinit':
        p.EFTCAMB.initialize_parameters(p, {'EFTflag': 1, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1,
                                            'EFTOmega0': 0.01, 'PureEFTmodelGamma1': 0, 'PureEFTmodelGamma2': 0,
                                            'PureEFTmodelGamma3': 0, 'PureEFTmodelGamma4': 0,
                                            'PureEFTmodelGamma5': 0, 'PureEFTmodelGamma6': 0,
                                            'EFTwDE': 0}, print_header=False)
        what = 'get_results'
show(p)
if what in ('get_results', 'get_background'):
    if kind not in ('valid',):
        try:
            p.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12)
        except TypeError:  # B1 (CAMB 1.6.5) + numpy 2: BBN YHe returns an array
            print('NOTE set_cosmology TypeError, retry with YHe=0.2454', flush=True)
            p.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12, YHe=0.2454)
    p.WantTransfer = False
    p.set_for_lmax(300, lens_potential_accuracy=0)
    try:
        if what == 'get_results':
            r = camb.get_results(p)
            cl = r.get_total_cls(300)
            print('VALUE get_results tau0 = %.6f  TT[200] = %.6e' % (r.tau0, cl[200, 0]), flush=True)
        else:
            r = camb.get_background(p)
            print('VALUE get_background H(z=1) = %.8f' % r.hubble_parameter(1.0), flush=True)
    except Exception as e:
        print('VALUE python exception', type(e).__name__, str(e)[:120], flush=True)
else:
    print('VALUE', what, '=', repr(getattr(p.EFTCAMB, what)())[:200], flush=True)
"""

GETTERS = ["model_name", "num_params", "param_names", "param_labels", "param_values"]
ALL_CASES = (
    [f"default_{g}" for g in GETTERS]
    + [f"setparamsgr_{g}" for g in GETTERS]
    + [f"valid_{g}" for g in GETTERS]
    + [f"failmodel_{g}" for g in GETTERS]
    + ["failmodel_get_background", "failmodel_get_results"]
    + [f"failomega_{g}" for g in ["num_params", "param_names", "param_values"]]
    + ["failomega_get_results", "valid_get_results", "default_get_results", "reinit_ok"]
)


def run_case(tree, case, cwd):
    env = dict(os.environ, OMP_NUM_THREADS="1", VERIFY_CWD=cwd)
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        r = subprocess.run(["nice", "-n", "10", sys.executable, "-c", CHILD, tree, case, INI], env=env,
                           capture_output=True, text=True, timeout=900)
    rc = r.returncode
    how = f"killed by {signal.Signals(-rc).name}" if rc < 0 else f"rc={rc}"
    keep = [ln for ln in r.stdout.splitlines() if ln.startswith(("VALUE", "STATE", "initialize", "NO EXC", "NOTE"))]
    err = [ln for ln in r.stderr.splitlines() if ln.strip()][-12:]
    return {"tree": tree, "case": case, "rc": rc, "how": how, "stdout": keep, "stderr_tail": err}


if __name__ == "__main__":
    out = sys.argv[1]
    args = sys.argv[2:]
    cases = ALL_CASES
    if "--cases" in args:
        i = args.index("--cases")
        cases = args[i + 1].split(",")
        args = args[:i] + args[i + 2:]
    cwd = os.path.dirname(os.path.abspath(out))
    os.makedirs(cwd, exist_ok=True)
    res = []
    for tree in args:
        for case in cases:
            d = run_case(tree, case, cwd)
            res.append(d)
            print(f"{os.path.basename(tree):6s} {case:26s} {d['how']:20s} {' | '.join(d['stdout'])}", flush=True)
            if d["rc"] != 0 and d["stderr_tail"]:
                for ln in d["stderr_tail"][-6:]:
                    print("        stderr:", ln[:160])
    with open(out, "w") as f:
        json.dump(res, f, indent=1)

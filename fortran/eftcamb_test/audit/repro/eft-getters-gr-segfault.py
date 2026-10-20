"""Repro (t9-wrapper): the EFTCAMB python getters num_params / param_names / param_labels / param_values kill the
interpreter (SIGSEGV) on any parameters object without an EFT model, e.g. a default camb.CAMBparams() (EFTflag = 0)
or one where initialize_parameters failed.
A failed initialize_parameters (ValueError) also leaves EFTflag /= 0 without a model (initialize_EFTCAMB_from_py,
11_EFTCAMB_python.f90:139-146, never resets it), and camb.get_results on that object then dies with SIGSEGV
(results.f90 CAMBdata_SetParams dereferences this%CP%EFTCAMB%model).

fortran/eftcamb/11_EFTCAMB_python.f90: EFTCAMB_get_num_params (l.200), EFTCAMB_get_param_names (l.214),
EFTCAMB_get_param_labels (l.228) and EFTCAMB_get_param_values (l.242) dereference EFTCAMB%model without checking
allocated(EFTCAMB%model) (EFTCAMB_get_model_name, l.255, does check). param_names/labels/values call num_params first.
Expected: 0 / [] for GR (as model_name() returns ''), or a python exception.

Run: /Users/marco/Software/miniforge3/bin/python eft-getters-gr-segfault.py [--tree PATH]
"""

import argparse
import os
import signal
import subprocess
import sys

TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")

CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
p = camb.CAMBparams()
if sys.argv[3] == 'failed_init':
    try:
        p.EFTCAMB.initialize_parameters(p, {'EFTflag': 1, 'PureEFTmodel': 9}, print_header=False)
    except ValueError as e:
        print('initialize_parameters raised', repr(e), flush=True)
print('EFTflag =', p.EFTCAMB.EFTflag, '| model_name() =', repr(p.EFTCAMB.model_name()), flush=True)
if sys.argv[2] == 'get_results':
    p.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12)
    p.set_for_lmax(300, lens_potential_accuracy=0)
    try:
        r = camb.get_results(p)
        print('value = results, tau0 =', r.tau0, flush=True)
    except Exception as e:
        print('value = python exception', repr(e)[:100], flush=True)
else:
    print('value =', getattr(p.EFTCAMB, sys.argv[2])(), flush=True)
"""

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=TREE)
    a = ap.parse_args()
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    for state in ("default", "failed_init"):
        for meth in ("model_name", "num_params", "param_names", "param_labels", "param_values", "get_results"):
            r = subprocess.run([sys.executable, "-c", CHILD, a.tree, meth, state], env=env, capture_output=True,
                               text=True)
            rc = r.returncode
            how = f"killed by {signal.Signals(-rc).name}" if rc < 0 else f"rc={rc}"
            last = [ln for ln in r.stdout.splitlines() if ln.startswith(("value", "EFTflag"))]
            print(f"{state:12s} {meth:13s}: {how:20s} {' | '.join(last)}   (expected rc=0 and a value)")

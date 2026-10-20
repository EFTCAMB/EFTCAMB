"""Independent check (verify) of finding initialize-parameters-bad-value-kills-python.

Bypasses camb.eftcamb.initialize_parameters: calls the two Fortran entry points of 11_EFTCAMB_python.f90 directly
through ctypes with hand-made key/value strings, so the Python-side str() conversion is not involved.
Each case runs in a child process; we record returncode, the forutils error line and whether the child reached
the line after the Fortran call.  Also runs controls (well-formed values), the same cases on B1 (pre-rebase
EFTCAMB) and an upstream-CAMB (B0) analogue (numpy/float values through camb.set_params).

Run: /Users/marco/Software/miniforge3/bin/python initparams_badvalue_indep.py
"""

import json
import os
import subprocess
import sys

TREES = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees"

# direct Fortran entry: key/value given verbatim as strings (no python conversion)
CHILD_DIRECT = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
from ctypes import c_int, byref, create_string_buffer
import camb
from camb import eftcamb as E
pairs = json.loads(sys.argv[2])
p = camb.CAMBparams()
p.set_cosmology(H0=67.5, YHe=0.245)  # explicit YHe: B1 BBN interpolation fails with numpy 2
E.EFTCAMB_initialize_read_params(create_string_buffer(b' '), c_int(1), create_string_buffer(b' '), c_int(1), c_int(0))
for k, v in pairs:
    kb, vb = k.encode(), v.encode()
    E.EFTCAMB_initialize_read_params(create_string_buffer(kb), c_int(len(kb)), create_string_buffer(vb), c_int(len(vb)),
                                     c_int(1))
err = c_int(0)
print('MARK before fortran init', flush=True)
E.EFTCAMB_initialize(byref(p), byref(err))
print('MARK returned eft_error=%d EFTflag=%d' % (err.value, p.EFTCAMB.EFTflag), flush=True)
"""

# upstream-CAMB analogue: numpy / float values for integer & logical fields through set_params
CHILD_B0 = r"""
import sys
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
kw = eval(sys.argv[2])
try:
    camb.set_params(H0=67.5, ombh2=0.022, omch2=0.12, **kw)
    print('MARK accepted', flush=True)
except Exception as e:
    print('MARK python exception', type(e).__name__, str(e)[:80], flush=True)
print('MARK alive', flush=True)
"""

BASE = [["EFTflag", "1"], ["PureEFTmodel", "1"], ["PureEFTmodelOmega", "1"], ["EFTOmega0", "0.01"], ["EFTwDE", "0"]]


def with_(key, val):
    d = dict(BASE)
    d[key] = val
    return [[k, v] for k, v in d.items()]


DIRECT_CASES = [
    ("control: all well-formed", BASE),
    ("control: logical 'T'", with_("EFT_ghost_stability", "T")),
    ("logical 'True'", with_("EFT_ghost_stability", "True")),
    ("integer '1.0' (EFTflag)", with_("EFTflag", "1.0")),
    ("integer '0.0' (EFTwDE)", with_("EFTwDE", "0.0")),
    ("double 'None'", with_("EFTOmega0", "None")),
    ("double \"b'0.01'\"", with_("EFTOmega0", "b'0.01'")),
    ("control: double '1e-2'", with_("EFTOmega0", "1e-2")),
]

B0_CASES = [
    "dict(lmax=np.float64(2500))",
    "dict(lmax=2500.0)",
    "dict(WantTensors=np.True_)",
    "dict(num_massive_neutrinos=np.float64(1))",
    "dict(halofit_version=None)",
]


def run(cmd, env):
    r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=300)
    out = r.stdout + r.stderr
    marks = [ln[5:] for ln in r.stdout.splitlines() if ln.startswith("MARK")]
    errs = [ln.strip() for ln in out.splitlines() if "Error for key" in ln or "STOP" in ln]
    return r.returncode, marks, errs


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import harness

    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    results = {}
    with harness.slot(1):
        # python-side conversion of the values in the finding
        import numpy as np

        conv = {
            "np.True_": str(np.True_), "1.0": str(1.0), "None": str(None), "b'0.01'": str(b"0.01"),
            "np.float64(0)": str(np.float64(0)), "np.int64(1)": str(np.int64(1)), "np.float64(0.01)": str(np.float64(0.01)),
        }
        print("python str() of values:", conv)
        results["conv"] = conv
        for tree in ("A0", "B1"):
            for label, pairs in DIRECT_CASES:
                rc, marks, errs = run([sys.executable, "-c", CHILD_DIRECT, f"{TREES}/{tree}", json.dumps(pairs)], env)
                print(f"[{tree} direct] {label:28s} rc={rc} marks={marks} errs={errs}")
                results[f"{tree}:{label}"] = dict(rc=rc, marks=marks, errs=errs)
        for kw in B0_CASES:
            rc, marks, errs = run([sys.executable, "-c", CHILD_B0, f"{TREES}/B0", kw], env)
            print(f"[B0 set_params] {kw:42s} rc={rc} marks={marks} errs={errs}")
            results[f"B0:{kw}"] = dict(rc=rc, marks=marks, errs=errs)
    out = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit/work/verify/initparams_badvalue_indep.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(results, f, indent=1)
    print("wrote", out)

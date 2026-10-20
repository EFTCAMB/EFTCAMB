"""Independent check: drive the ctypes entry points directly (bypassing initialize_parameters' str() conversion)
with raw strings, on A0 and B1; also check controls (np.int64, good values) and whether interpreter survives."""
import os, subprocess, sys
BASE = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/"
CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
from camb import eftcamb as E
from ctypes import c_int, byref, create_string_buffer
raw = dict(EFTflag='1', PureEFTmodel='1', PureEFTmodelOmega='1', EFTOmega0='0.01', EFTwDE='0')
raw[sys.argv[2]] = sys.argv[3]
p = camb.CAMBparams(); p.set_cosmology(H0=67.5)
f = E.EFTCAMB_initialize_read_params
f(create_string_buffer(b' '), c_int(1), create_string_buffer(b' '), c_int(1), c_int(0))
for k, v in raw.items():
    k, v = k.encode(), v.encode()
    f(create_string_buffer(k), c_int(len(k)), create_string_buffer(v), c_int(len(v)), c_int(1))
t = c_int(0)
try:
    E.EFTCAMB_initialize(byref(p), byref(t))
    print('RESULT returned err=%d model=%s' % (t.value, p.EFTCAMB.EFTflag), flush=True)
except Exception as e:
    print('RESULT exc', e, flush=True)
"""
CASES = [("EFTflag", "1"), ("EFTflag", "1.0"), ("EFTwDE", "0.0"), ("EFTOmega0", "None"), ("EFT_ghost_stability", "True"),
         ("EFT_ghost_stability", "T")]
env = dict(os.environ, OMP_NUM_THREADS="1"); env.pop("PYTHONPATH", None)
for tree in ("A0", "B1"):
    for k, v in CASES:
        r = subprocess.run([sys.executable, "-c", CHILD, BASE + tree, k, v], env=env, capture_output=True, text=True)
        out = r.stdout + r.stderr
        res = [l for l in out.splitlines() if l.startswith("RESULT") or "Error for key" in l or "STOP" in l]
        print(f"{tree} {k}={v!r:8s} rc={r.returncode} {res}")
# controls through python API on A0
CTRL = r"""
import sys; sys.path.insert(0, sys.argv[1]); import numpy as np, camb
d = {'EFTflag': np.int64(1), 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1, 'EFTOmega0': np.float64(0.01), 'EFTwDE': 0,
     'EFT_ghost_stability': True}
p = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.12, **d); print('RESULT ok', p.EFTCAMB.EFTflag, flush=True)
"""
r = subprocess.run([sys.executable, "-c", CTRL, BASE + "A0"], env=env, capture_output=True, text=True)
print("control np.int64/np.float64/bool:", r.returncode, [l for l in r.stdout.splitlines() if "RESULT" in l], r.stderr[-300:])

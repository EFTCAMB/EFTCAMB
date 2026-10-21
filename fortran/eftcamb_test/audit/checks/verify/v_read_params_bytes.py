"""Independent check of read-parameters-roundtrip-bytes: (1) python semantics of float/int on bytes,
(2) initialize_parameters with str values '1.d-4' and plain output_root (bypassing read_parameters) works,
(3) decoding read_parameters output with a d->e conversion restores a faithful round trip."""
import os, subprocess, sys
REPO = "/Users/marco/IW/eftcamb_update/EFTCAMB"
TREE = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0"
for s in (b"1.d-4", b"1e-4", b"5"):
    for fn in (int, float):
        try: print(f"{fn.__name__}({s!r}) = {fn(s)!r}")
        except ValueError as e: print(f"{fn.__name__}({s!r}) -> ValueError")
CHILD = r"""
import sys, re
sys.path.insert(0, sys.argv[1]); import camb
mode = sys.argv[3]
p = camb.read_ini(sys.argv[2]); d = p.EFTCAMB.read_parameters()
if mode == 'manual':
    d = {k: v for k, v in d.items()}
    d['EFTB0'] = '1.d-4'; d['output_root'] = 'abc'
elif mode == 'decoded':
    d = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in d.items()}
q = camb.CAMBparams(); q.set_cosmology(H0=67)
q.EFTCAMB.initialize_parameters(q, d, print_header=False)
d2 = q.EFTCAMB.read_parameters()
print('EFTB0 after:', repr(d2.get('EFTB0')), 'output_root after:', repr(d2.get('output_root')))
if mode == 'decoded':
    d0 = p.EFTCAMB.read_parameters()
    diff = [k for k in d0 if d0[k] != d2.get(k)]
    print('keys differing after decoded round trip:', diff, 'nkeys', len(d0), len(d2))
"""
env = dict(os.environ, OMP_NUM_THREADS="1"); env.pop("PYTHONPATH", None)
for model in ("4_Des_F_R_LCDM_2.ini", "2_PEFT_Omega_power_1.ini"):
    for mode in ("manual", "decoded"):
        r = subprocess.run([sys.executable, "-c", CHILD, TREE, f"{REPO}/fortran/eftcamb_test/parameters/{model}", mode],
                           env=env, capture_output=True, text=True, cwd=f"{REPO}/fortran")
        out = [l for l in (r.stdout + r.stderr).splitlines() if l.startswith(("EFTB0", "keys", "Error", "ERROR", "Trace"))]
        print(model, mode, "rc", r.returncode, out)

"""Repro (track 6): the stability cache carries 3x the massive-neutrino pressure.

EFTStabilityComputation (fortran/eftcamb/09_EFTCAMB_stability.f90) fills the massive-neutrino totals three times.
The first loop (l.461-471) resets gpinu_tot, but the second (l.488-512) and third (l.517-541) loops reset
gpinudot_tot twice and never reset gpinu_tot before accumulating `gpinu_tot = gpinu_tot + grhormass_t*gpinu`.
So compute_H_derivs (called after loop 2 and after loop 3) sees 2x and then 3x the massive-neutrino pressure.
Models whose Hdotdot / Hdotdotdot use eft_cache%gpinu_tot (Pure EFT and the mc-quintessence designer that extends it,
Horava, designer 3+1) then feed wrong Hdotdot / Hdotdotdot into the pi-equation coefficients (mass_math check), the
mass eigenvalues (mass check) and the positivity / Minkowski coefficients.  The ghost and gradient coefficients do
not use Hdotdot.

Reference: the python bridge get_eft_functions (11_EFTCAMB_python.f90) fills the same cache with a single reset loop.
Both are evaluated here for the same CAMBdata object at the same scale factors.

Run: /Users/marco/Software/miniforge3/bin/python stability-gpinu-tot-tripled.py [tree [fixed_tree]]
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def _run(*args, **kw):
    """subprocess.run holding one machine-wide audit CPU slot (if the audit harness is available), at nice 10."""
    import contextlib
    import os as _os
    import sys as _sys
    from pathlib import Path as _P
    _sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
    try:
        import harness as _h
        ctx = _h.slot(1)
    except Exception:  # noqa: BLE001
        ctx = contextlib.nullcontext()
    with ctx:
        return subprocess.run(*args, preexec_fn=lambda: _os.nice(10), **kw)


ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREE = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/trees/A0"
PAR = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")
MODEL = "2_PEFT_gamma3_const_1"

WORKER = r"""
import ctypes, json, sys
from ctypes import byref, c_double, c_int
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
from camb.baseconfig import camblib as lib
from camb.eftcamb import EFTCAMB_timestep_cache
cp = camb.read_ini(sys.argv[2], no_validate=True)
cp.EFTCAMB.EFTCAMB_skip_stability = True
res = camb.get_background(cp)
P = res.Params
a = np.array([1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1.0])
fields, br = cp.EFTCAMB.get_eft_functions(res, a)
test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
buf = (ctypes.c_char * (1 << 18))()
names = ["gpinu_tot", "Hdotdot", "Hdotdotdot", "EFTpiC", "posbound2", "tachpos", "EFT_kinetic", "EFT_gradient"]
off = {n: getattr(EFTCAMB_timestep_cache, n).offset for n in names}
out = {"camb": camb.__file__, "a": a.tolist(), "num_nu_massive": int(P.num_nu_massive), "omnuh2": float(P.omnuh2),
       "stab": {n: [] for n in names}, "bridge": {n: br[n].tolist() for n in names}}
for x in a:
    test(byref(c_double(x)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf)
    raw = bytes(buf[: max(off.values()) + 8])
    for n in names:
        out["stab"][n].append(float(np.frombuffer(raw, dtype=np.float64, count=1, offset=off[n])[0]))
print("RESULT " + json.dumps(out))
"""

def run(tree):
    with tempfile.TemporaryDirectory() as d:
        ini = Path(d) / "m.ini"
        ini.write_text(f"DEFAULT({PAR / (MODEL + '.ini')})\noutput_root = {d}/out\nfeedback_level = 0\n"
                       "EFTCAMB_write_background = F\n")
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        p = _run([sys.executable, "-c", WORKER, tree, str(ini)], capture_output=True, text=True, cwd=d, env=env)
        res = [json.loads(line[7:]) for line in p.stdout.splitlines() if line.startswith("RESULT ")]
        if not res:
            print(p.stdout[-2000:], p.stderr[-2000:])
            sys.exit(1)
        return res[0]


trees = [TREE] + sys.argv[2:3]
runs = [run(t) for t in trees]
for t, r in zip(trees, runs):
    print("tree:", t)
    print(f"  model {MODEL}: num_nu_massive = {r['num_nu_massive']}, omnuh2 = {r['omnuh2']}")
    print("       a    | gpinu_tot stability/bridge | Hdotdot (stab-bridge)/|bridge| | Hdotdotdot (stab-bridge)/|bridge|")
    s, b = r["stab"], r["bridge"]
    for i, a in enumerate(r["a"]):
        print(f"    {a:8.1e} |          {s['gpinu_tot'][i] / b['gpinu_tot'][i]:7.4f}           |"
              f"           {(s['Hdotdot'][i] - b['Hdotdot'][i]) / abs(b['Hdotdot'][i]): 10.2e}           |"
              f"           {(s['Hdotdotdot'][i] - b['Hdotdotdot'][i]) / abs(b['Hdotdotdot'][i]): 10.2e}")
if len(runs) == 2:
    print("downstream stability coefficients, relative change first tree -> second tree:")
    print("       a    |   EFTpiC    posbound2   tachpos   EFT_kinetic EFT_gradient")
    s0, s1 = runs[0]["stab"], runs[1]["stab"]
    for i, a in enumerate(runs[0]["a"]):
        print(f"    {a:8.1e} |" + "".join(
            f" {(s0[n][i] - s1[n][i]) / abs(s1[n][i]) if s1[n][i] else float('nan'): 10.2e}"
            for n in ("EFTpiC", "posbound2", "tachpos", "EFT_kinetic", "EFT_gradient")))
print("EXPECTED: gpinu_tot ratio = 1 and Hdotdot/Hdotdotdot equal to the bridge values (same model, a, CAMBdata).")

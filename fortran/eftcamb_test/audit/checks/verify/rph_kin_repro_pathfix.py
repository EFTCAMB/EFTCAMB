"""Repro (track 6): spurious scalar-ghost rejection of a healthy RPH model at early times from catastrophic
cancellation in EFT_kinetic (06_abstract_EFTCAMB_model.f90:690).

Model: the test ini 3_RPH_MB_OmegaDE_1 (alpha_M = Omega_DE, alpha_B = -Omega_DE, alpha_K = alpha_T = 0, LCDM
background).  Analytic no-ghost (EFTCAMB alpha_B convention): D = alpha_K + 6 alpha_B^2 = 6 Omega_DE^2 > 0 at all times,
EFT_kinetic = 18 M*^6/m0^6 H^2 D > 0.  RPH builds gamma1 from alpha_K and c (4 gamma1 H0^2 a^2 = alpha_K M H^2 - 2c),
so EFT_kinetic = 9 M (4 c M + 8 gamma1 M a^2 H0^2 + ...) cancels 4cM - 4cM; at a < 1e-5, |c| exceeds the true kinetic
term by > 1e16 and its sign is roundoff.  The test ini hides this by setting EFTCAMB_stability_time = 0.1; with the
code default (1e-10) the model is rejected.

Run: /Users/marco/Software/miniforge3/bin/python rph-kinetic-early-roundoff-ghost.py [tree]
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


ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREE = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/trees/A0"
PAR = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"

WORKER = r"""
import sys, json, ctypes
from ctypes import byref, c_double, c_int
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
from camb.baseconfig import camblib as lib
from camb.eftcamb import EFTCAMB_timestep_cache as C

out = {"verdict": {}}
for st in (0.1, 1e-10):
    cp = camb.read_ini(sys.argv[2], no_validate=True)
    cp.EFTCAMB.EFTCAMB_stability_time = st
    try:
        camb.get_background(cp); out["verdict"][str(st)] = "stable"
    except Exception as e:
        out["verdict"][str(st)] = "unstable"
cp = camb.read_ini(sys.argv[2], no_validate=True)
cp.EFTCAMB.EFTCAMB_skip_stability = True
P = camb.get_background(cp).Params
test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
buf = (ctypes.c_char * (1 << 18))()
grid = np.logspace(-10, 0, 2001)
kin, ana, cterm = [], [], []
for a in grid:
    test(byref(c_double(a)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf)
    raw = bytes(buf[:8 * 200])
    g = lambda f: float(np.frombuffer(raw, np.float64, 1, getattr(C, f).offset)[0])
    M = g("Meff2")
    kin.append(g("EFT_kinetic"))
    ana.append(18 * M**3 * g("adotoa")**2 * (g("alphaK") + 6 * g("alphaB")**2))
    cterm.append(36 * M**2 * g("EFTc"))
kin, ana, cterm = map(np.array, (kin, ana, cterm))
neg = kin < 0
out["n_points"] = len(grid)
out["n_code_kinetic_negative"] = int(neg.sum())
out["n_analytic_negative"] = int((ana < 0).sum())
out["a_range_code_negative"] = [float(grid[neg].min()), float(grid[neg].max())] if neg.any() else None
out["rows"] = [{"a": float(grid[i]), "kin_code": float(kin[i]), "kin_analytic": float(ana[i]), "36M^2c": float(cterm[i])}
               for i in (100, 300, 500, 700, 1000, 1400, 2000)]
print("RESULT " + json.dumps(out))
"""

with tempfile.TemporaryDirectory() as d:
    ini = Path(d) / "m.ini"
    ini.write_text(f"DEFAULT({PAR}/3_RPH_MB_OmegaDE_1.ini)\nfeedback_level = 0\nEFTCAMB_write_background = F\n"
                   f"output_root = {d}/out\n")
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    p = _run([sys.executable, "-c", WORKER, TREE, str(ini)], capture_output=True, text=True, cwd=d, env=env)
    res = [json.loads(l[7:]) for l in p.stdout.splitlines() if l.startswith("RESULT ")]
    if not res:
        print(p.stdout[-2000:], p.stderr[-2000:])
        sys.exit(1)
    r = res[0]
print("tree:", TREE)
print("verdict with EFTCAMB_stability_time = 0.1 (test ini):", r["verdict"]["0.1"])
print("verdict with EFTCAMB_stability_time = 1e-10 (code default):", r["verdict"]["1e-10"])
print(f"EFT_kinetic < 0 at {r['n_code_kinetic_negative']}/{r['n_points']} log-spaced points in [1e-10, 1], "
      f"a in {r['a_range_code_negative']}; analytic 18 M^3 H^2 (alpha_K + 6 alpha_B^2) < 0 at {r['n_analytic_negative']}")
for x in r["rows"]:
    print(f"  a={x['a']:.1e} EFT_kinetic code={x['kin_code']: .3e} analytic={x['kin_analytic']: .3e} "
          f"size of cancelling term 36 M^2 EFTc={x['36M^2c']: .3e}")
print("EXPECTED: stable for both stability times (D = 6 Omega_DE^2 > 0).")

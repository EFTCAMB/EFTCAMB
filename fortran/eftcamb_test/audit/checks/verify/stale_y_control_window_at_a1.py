"""Repro (track 6): the third sampling loop of EFTCAMB_Stability_Check ("log sampling close to aend") never updates
y, so it evaluates a = astart 1000 times.  An instability confined to a narrow interval just below a = 1 is missed.

Model: Pure EFT, LCDM background, gamma1(a) = -1e-3 + 2.5e4 (a - 0.9995)^2 (Taylor order 2).  The scalar kinetic term
is 72 a^2 H0^2 gamma1, negative exactly for |a - 0.9995| < 2e-4.  Expected: theory unstable.

Run: /Users/marco/Software/miniforge3/bin/python stability-sampler-stale-y.py [tree]
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
    _sys.path.insert(0, str(_P(__file__).resolve().parents[1] / "checks"))
    try:
        import harness as _h
        ctx = _h.slot(1)
    except Exception:  # noqa: BLE001
        ctx = contextlib.nullcontext()
    with ctx:
        return subprocess.run(*args, preexec_fn=lambda: _os.nice(10), **kw)


ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREE = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/trees/A0"
BASE = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/base_params.ini"

INI = f"""
EFTflag = 1
PureEFTmodel = 1
PureEFTHorndeski = F
EFTwDE = 0
PureEFTmodelOmega = 0
PureEFTmodelGamma1 = 5
PureEFTmodelGamma2 = 0
PureEFTmodelGamma3 = 0
PureEFTmodelGamma4 = 0
PureEFTmodelGamma5 = 0
PureEFTmodelGamma6 = 0
EFTGamma1_Taylor_order = 2
EFTGamma1a0 = 0.99995
EFTGamma10 = -1.e-3
EFTGamma11 = 0.
EFTGamma12 = 2.5e4
EFT_ghost_stability = T
EFT_gradient_stability = T
EFT_mass_stability = F
EFTCAMB_stability_time = 1.e-10
feedback_level = 0
DEFAULT({BASE})
"""

WORKER = r"""
import sys, json, ctypes
from ctypes import byref, c_double, c_int
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
from camb.baseconfig import camblib as lib
from camb.eftcamb import EFTCAMB_timestep_cache

ini = sys.argv[2]
out = {"camb": camb.__file__}
# 1) full code path, stability on
cp = camb.read_ini(ini, no_validate=True)
try:
    camb.get_background(cp)
    out["code_verdict"] = "stable"
except Exception as e:
    out["code_verdict"] = "unstable: " + str(e).splitlines()[-1]
# 2) per-point EFTTestStability on a fine grid (stability skipped in the background init)
cp = camb.read_ini(ini, no_validate=True)
cp.EFTCAMB.EFTCAMB_skip_stability = True
res = camb.get_background(cp)
P = res.Params
test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
check = getattr(lib, "__eftcamb_stability_MOD_eftcamb_stability_check")
buf = (ctypes.c_char * (1 << 18))()
off = EFTCAMB_timestep_cache.EFT_kinetic.offset
grid = np.linspace(0.999, 1.0, 1001)
bad = []
for a in grid:
    ok = test(byref(c_double(a)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf)
    kin = np.frombuffer(bytes(buf[off:off + 8]), dtype=np.float64)[0]
    if not ok:
        bad.append((float(a), float(kin)))
out["fine_grid_unstable_points"] = len(bad)
out["fine_grid_unstable_range"] = [bad[0][0], bad[-1][0]] if bad else None
out["kinetic_at_0.9995"] = [k for a, k in bad if abs(a - 0.9995) < 1e-9]
# 3) the sampler itself on the same objects
succ = c_int(0)
check(byref(succ), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), byref(c_double(1e-10)), byref(c_double(1.0)),
      byref(c_double(10.0)))
out["sampler_verdict"] = "stable" if succ.value else "unstable"
# 4) the a values the three loops visit (formulas copied from 09_EFTCAMB_stability.f90)
i = np.arange(1000)
lin = 1e-10 + i * (1 - 1e-10) / 999
y = -10 + i * 10 / 999
log1 = 1e-10 + (1 - 1e-10) * 10**y
log3_code = np.full(1000, 1 + (1e-10 - 1) * 10**y[-1])   # y is never updated in loop 3: stays 0
log3_intended = 1 + (1e-10 - 1) * 10**y
win = lambda g: int(np.sum((g > 0.9993) & (g < 0.9997)))
out["points_in_unstable_window"] = {"loop1_linear": win(lin), "loop2_log_astart": win(log1),
                                    "loop3_as_coded": win(log3_code), "loop3_intended": win(log3_intended),
                                    "loop3_as_coded_unique_a": sorted(set(log3_code.tolist()))}
print("RESULT " + json.dumps(out))
"""

with tempfile.TemporaryDirectory() as d:
    ini = Path(d) / "window.ini"
    ini.write_text(INI)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    p = _run([sys.executable, "-c", WORKER, TREE, str(ini)], capture_output=True, text=True, cwd=d, env=env)
    res = [json.loads(l[7:]) for l in p.stdout.splitlines() if l.startswith("RESULT ")]
    if not res:
        print(p.stdout[-2000:], p.stderr[-2000:])
        sys.exit(1)
    r = res[0]
print("tree:", TREE)
print("per-point EFTTestStability on 1001 points in [0.999, 1]: unstable points =", r["fine_grid_unstable_points"],
      "range =", r["fine_grid_unstable_range"], " (analytic window: 0.9993 < a < 0.9997)")
print("EFT_kinetic at a = 0.9995:", r["kinetic_at_0.9995"], " (analytic 72 H0^2 a^2 gamma1 < 0)")
print("sampler grid points inside the window:", json.dumps(r["points_in_unstable_window"]))
print("OBSERVED: EFTCAMB_Stability_Check ->", r["sampler_verdict"], "; get_background ->", r["code_verdict"])
print("EXPECTED: unstable (ghost for 0.9993 < a < 0.9997; the intended loop 3 has",
      r["points_in_unstable_window"]["loop3_intended"], "points there)")

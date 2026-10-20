"""Independent consequence check for finding stability-third-pass-stale-y.

A Pure-EFT model on an LCDM background with gamma1(a) = g0 + g2 (a - a0)^2 (Taylor, order 2), g0 < 0 < g2, has a
negative scalar kinetic term (ghost) only for |a - a0| < sqrt(-g0/g2).  We place that window just below aend, in a gap
of the linear pass (spacing ~(aend-astart)/999) and of the log-near-astart pass, so only the intended third pass
(log sampling close to aend) can see it.

For each tree (A0 = code under audit, FIXY = A0 + one-line y update in loop 3) we record:
  1. full code path: camb.get_background with the stability check on (aend = 1), feedback_level 3 (captures
     "Instability detected at a =" lines);
  2. direct ctypes call of EFTCAMB_Stability_Check on the initialised objects with a chosen [astart, aend];
  3. sanity: EFTTestStability evaluated directly at a = a0 (must be unstable) and at aend (must be stable).

Usage: python stab3_window_blackbox.py
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45"
             "/scratchpad/audit")
REPO = Path(__file__).resolve().parents[5]
BASE = REPO / "fortran/eftcamb_test/parameters/base_params.ini"
TREES = {"A0": AUDIT / "trees/A0", "FIXY": AUDIT / "work/t6-stability/fixy"}
OUT = AUDIT / "work/verify/stab3"

# cases: (name, a0, half-width, g0, astart, aend for the direct call)
CASES = [
    ("win_0.99985_aend1", 0.99985, 5e-5, -1e-4, 1e-10, 1.0),
    ("win_0.49985_aend0.5", 0.49985, 5e-5, -1e-4, 1e-10, 0.5),
    ("win_0.99995_astart0.9", 0.99995, 2e-5, -1e-4, 0.9, 1.0),
]


def ini_text(a0, hw, g0):
    g2 = -g0 / hw**2
    return f"""EFTflag = 1
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
EFTGamma1a0 = {a0!r}
EFTGamma10 = {g0!r}
EFTGamma11 = 0.
EFTGamma12 = {g2!r}
EFT_ghost_stability = T
EFT_gradient_stability = T
EFT_mass_stability = F
EFTCAMB_stability_time = 1.e-10
feedback_level = 3
EFTCAMB_write_background = F
DEFAULT({BASE})
"""


WORKER = r"""
import sys, json, ctypes
from ctypes import byref, c_double, c_int
sys.path.insert(0, sys.argv[1])
import camb
from camb.baseconfig import camblib as lib
ini, a0, astart, aend = sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
out = {"camb": camb.__file__}
cp = camb.read_ini(ini, no_validate=True)
try:
    camb.get_background(cp)
    out["full_path"] = "stable"
except Exception as e:
    out["full_path"] = "unstable: " + str(e).strip().splitlines()[-1]
sys.stdout.flush()
cp = camb.read_ini(ini, no_validate=True)
cp.EFTCAMB.EFTCAMB_skip_stability = True
cp.EFTCAMB.EFTCAMB_feedback_level = 0
res = camb.get_background(cp)
P = res.Params
check = lib.__eftcamb_stability_MOD_eftcamb_stability_check
test = lib.__eftcamb_stability_MOD_eftteststability
test.restype = c_int
buf = (ctypes.c_char * (1 << 18))()
out["direct_test_at_a0"] = bool(test(byref(c_double(a0)), byref(c_double(10.0)), byref(P.EFTCAMB),
                                     byref(P.EFTCAMB_parameter_cache), buf))
out["direct_test_at_aend"] = bool(test(byref(c_double(aend)), byref(c_double(10.0)), byref(P.EFTCAMB),
                                       byref(P.EFTCAMB_parameter_cache), buf))
succ = c_int(-7)
check(byref(succ), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), byref(c_double(astart)),
      byref(c_double(aend)), byref(c_double(10.0)))
out["direct_check"] = "stable" if succ.value else "unstable"
print("RESULT " + json.dumps(out), flush=True)
"""


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    results = {}
    with tempfile.TemporaryDirectory(dir=OUT) as d:
        for name, a0, hw, g0, astart, aend in CASES:
            ini = Path(d) / f"{name}.ini"
            ini.write_text(ini_text(a0, hw, g0))
            for tname, tree in TREES.items():
                env = dict(os.environ, OMP_NUM_THREADS="1")
                env.pop("PYTHONPATH", None)
                with harness.slot(1):
                    p = subprocess.run([sys.executable, "-c", WORKER, str(tree), str(ini), str(a0), str(astart),
                                        str(aend)], capture_output=True, text=True, cwd=d, env=env, timeout=500,
                                       preexec_fn=lambda: os.nice(10))
                r = [json.loads(line[7:]) for line in p.stdout.splitlines() if line.startswith("RESULT ")]
                det = [line.strip() for line in p.stdout.splitlines() if "Instability detected" in line]
                why = [line.strip() for line in p.stdout.splitlines() if "Physical instability" in line]
                rec = r[0] if r else {"error": p.stdout[-1500:] + p.stderr[-1500:]}
                rec["detected_lines"] = det
                rec["reason_lines"] = why[:3]
                results[f"{name}|{tname}"] = rec
                print(f"{name:26s} {tname:5s} full_path(aend=1)={rec.get('full_path')!s:40s} "
                      f"direct_check[{astart:g},{aend:g}]={rec.get('direct_check')}  "
                      f"test(a0)={rec.get('direct_test_at_a0')} test(aend)={rec.get('direct_test_at_aend')}  "
                      f"{det[:1]}", flush=True)
    (OUT / "window_blackbox.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()

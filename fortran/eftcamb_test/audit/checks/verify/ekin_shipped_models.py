"""For every shipped model, evaluate the stability-path timestep cache (EFTTestStability, which runs
EFTStabilityComputation: background, pi factors, stability factors) at a = 0.1 ... 1 on the A0 tree and compare the
code's EFT_kinetic with
  (i)  the documented formula of the numerical notes (single h0 on the cross term == proposed fix), evaluated from the
       code's own EFT functions, and
  (ii) 36 Mt^2 EFTpiA1, EFTpiA1 being the kinetic coefficient the code uses in the pi equation
       (identity valid for gamma3 = -gamma4, gamma6 = 0, where the ghost check is applied).
(The get_eft_functions output in the reference runs has EFT_kinetic = EFTpiA1 = 0, so it cannot be used.)
Each model runs in its own subprocess holding one harness CPU slot.
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ROOT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45"
            "/scratchpad/audit")
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "trees/A0"
OUT = ROOT / f"work/verify/ekin/shipped_models_{TREE.name}.json"

WORKER = r"""
import sys, json, ctypes
from ctypes import byref, c_double, c_int
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
from camb.baseconfig import camblib as lib
from camb.eftcamb import EFTCAMB_timestep_cache as C
cp = camb.read_ini(sys.argv[2], no_validate=True)
try:  # the stability module's saved state is only set up by a run with the stability check (as in the t6 repro)
    camb.get_background(cp); verdict = "stable"
except Exception:
    verdict = "unstable"
cp = camb.read_ini(sys.argv[2], no_validate=True)
cp.EFTCAMB.EFTCAMB_skip_stability = True
P = camb.get_background(cp).Params
test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
buf = (ctypes.c_char * (1 << 18))()
h0 = P.EFTCAMB_parameter_cache.h0_Mpc
rows = []
for a in np.linspace(0.1, 1.0, 10):
    test(byref(c_double(a)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf)
    raw = bytes(buf[:8 * 400])
    g = lambda f: float(np.frombuffer(raw, np.float64, 1, getattr(C, f).offset)[0])
    H, OP, Om, c = g("adotoa"), g("EFTOmegaP"), g("EFTOmegaV"), g("EFTc")
    g1, g2, g3, g4, g6 = g("EFTGamma1V"), g("EFTGamma2V"), g("EFTGamma3V"), g("EFTGamma4V"), g("EFTGamma6V")
    Mt = 1 + Om - g4
    fixed = 9*Mt*(4*c*Mt + 3*H**2*OP**2*a**2 + a**2*h0*(h0*(3*g2**2 + 8*g1*Mt) + 6*H*g2*OP))
    rows.append({"a": float(a), "code": g("EFT_kinetic"), "fixed": fixed, "A1x36Mt2": 36*Mt**2*g("EFTpiA1"),
                 "glpv": (g3 + g4 == 0) and g6 == 0, "cross": 54*Mt*a**2*h0*H*g2*OP})
print("RESULT " + json.dumps({"h0": h0, "verdict": verdict, "rows": rows}))
"""

res = {}
for ini in harness.model_inis():
    with tempfile.TemporaryDirectory() as d:
        w = Path(d) / "m.ini"
        w.write_text(f"DEFAULT({ini})\nfeedback_level = 0\nEFTCAMB_write_background = F\n")
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        with harness.slot(1):
            p = subprocess.run([sys.executable, "-c", WORKER, str(TREE), str(w)], capture_output=True, text=True,
                               cwd=d, env=env, timeout=600, preexec_fn=lambda: os.nice(10))
    lines = [ln for ln in p.stdout.splitlines() if ln.startswith("RESULT ")]
    if not lines:
        res[ini.stem] = {"error": (p.stdout + p.stderr)[-300:]}
        print(f"{ini.stem:34s} ERROR rc={p.returncode}")
        continue
    r = json.loads(lines[0][7:])
    rows = [x for x in r["rows"] if x["fixed"] != 0]
    rel = max((abs(x["code"] / x["fixed"] - 1) for x in rows), default=0.0)
    flips = sum((x["code"] > 0) != (x["fixed"] > 0) for x in rows)
    glpv = all(x["glpv"] for x in r["rows"])
    a1 = [x["A1x36Mt2"] / x["fixed"] for x in rows] if glpv else []
    a1c = [x["A1x36Mt2"] / x["code"] for x in rows if x["code"] != 0] if glpv else []
    res[ini.stem] = {"rel_code_vs_fixed": rel, "sign_flips": flips, "glpv": glpv, **r}
    print(f"{ini.stem:34s} {r['verdict']:8s} glpv={glpv!s:5s} max|code/fixed-1|={rel:.3e} sign_flips={flips}"
          + (f" 36Mt2A1/fixed=[{min(a1):.6f},{max(a1):.6f}] 36Mt2A1/code=[{min(a1c):.3f},{max(a1c):.3f}]" if a1 else ""),
          flush=True)
OUT.write_text(json.dumps(res, indent=1))
bad = [k for k, v in res.items() if v.get("rel_code_vs_fixed", 0) > 1e-6]
print(f"models where code != documented formula (>1e-6): {len(bad)} / {len(res)}; "
      f"with a sign flip on the grid: {[k for k, v in res.items() if v.get('sign_flips')]}")

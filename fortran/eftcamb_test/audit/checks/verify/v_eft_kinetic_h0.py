"""Repro (track 6): EFT_kinetic (scalar no-ghost coefficient, 06_abstract_EFTCAMB_model.f90:690) multiplies the cross
term 6 adotoa gamma2 Omega' by h0_Mpc^2 instead of h0_Mpc (dimensionally inconsistent: Mpc^-3 vs Mpc^-2).

Independent check: for RPH (alpha basis, EFTCAMB convention) the EFT -> alpha mapping gives exactly
    EFT_kinetic = 18 M*^6/m0^6 H^2 (alpha_K + 6 alpha_B^2)            (H = conformal Hubble, M*^2/m0^2 = Meff2)
which the code reproduces when Omega' = 0, but not with a running Planck mass.

Model A: the test ini 3_RPH_MK_OmegaDE_1 (alpha_M = 0.5 Omega_DE, alpha_K = 0.1 Omega_DE, alpha_B = 0).
Model B: same alpha_M, alpha_K = -0.2 Omega_DE^2 (a scalar ghost at every time: D = alpha_K < 0).

Run: /Users/marco/Software/miniforge3/bin/python eft-kinetic-cross-term-extra-h0.py [tree]
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

INI_A = f"DEFAULT({PAR}/3_RPH_MK_OmegaDE_1.ini)\nfeedback_level = 0\nEFTCAMB_write_background = F\n"
INI_B = INI_A + "RPHkineticitymodel_ODE = 3\nRPHkineticity_ODE0 = -0.2\nRPHkineticity_ODEExp = 2.0\n"

WORKER = r"""
import sys, json, ctypes
from ctypes import byref, c_double, c_int
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
from camb.baseconfig import camblib as lib
from camb.eftcamb import EFTCAMB_timestep_cache as C

out = {}
for tag, ini in (("A", sys.argv[2]), ("B", sys.argv[3])):
    cp = camb.read_ini(ini, no_validate=True)
    try:
        camb.get_background(cp); verdict = "stable"
    except Exception as e:
        verdict = "unstable"
    cp = camb.read_ini(ini, no_validate=True)
    cp.EFTCAMB.EFTCAMB_skip_stability = True
    P = camb.get_background(cp).Params
    test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
    buf = (ctypes.c_char * (1 << 18))()
    rows = []
    for a in (0.1, 0.3, 0.5, 0.7, 1.0):
        test(byref(c_double(a)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf)
        raw = bytes(buf[:8 * 200])
        g = lambda f: float(np.frombuffer(raw, np.float64, 1, getattr(C, f).offset)[0])
        H, M, aK, aB = g("adotoa"), g("Meff2"), g("alphaK"), g("alphaB")
        h0 = P.EFTCAMB_parameter_cache.h0_Mpc
        Mt = 1 + g("EFTOmegaV") - g("EFTGamma4V")
        c, g1, g2, OP = g("EFTc"), g("EFTGamma1V"), g("EFTGamma2V"), g("EFTOmegaP")
        fixed = 9*Mt*(4*c*Mt + 3*H**2*OP**2*a**2 + a**2*h0*(h0*(3*g2**2 + 8*g1*Mt) + 6*H*g2*OP))
        rows.append({"a": a, "alpha_K": aK, "alpha_M": g("alphaM"), "D=alpha_K+6alpha_B^2": aK + 6*aB**2,
                     "EFT_kinetic_code": g("EFT_kinetic"), "analytic_18M^3H^2D": 18*M**3*H**2*(aK + 6*aB**2),
                     "formula_with_single_h0": fixed})
    out[tag] = {"verdict": verdict, "stability_time": float(P.EFTCAMB.EFTCAMB_stability_time), "rows": rows}
print("RESULT " + json.dumps(out))
"""

with tempfile.TemporaryDirectory() as d:
    a, b = Path(d) / "A.ini", Path(d) / "B.ini"
    a.write_text(INI_A)
    b.write_text(INI_B)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    p = _run([sys.executable, "-c", WORKER, TREE, str(a), str(b)], capture_output=True, text=True, cwd=d,
                       env=env)
    res = [json.loads(l[7:]) for l in p.stdout.splitlines() if l.startswith("RESULT ")]
    if not res:
        print(p.stdout[-2000:], p.stderr[-2000:])
        sys.exit(1)
    r = res[0]
print("tree:", TREE)
for tag, label in (("A", "3_RPH_MK_OmegaDE_1 (alpha_K = 0.1 Omega_DE)"), ("B", "alpha_K = -0.2 Omega_DE^2 (ghost)")):
    print(f"Model {tag}: {label}; stability_time = {r[tag]['stability_time']}")
    for row in r[tag]["rows"]:
        print(f"  a={row['a']:.1f} alpha_K={row['alpha_K']: .4f} alpha_M={row['alpha_M']:.4f} "
              f"EFT_kinetic code={row['EFT_kinetic_code']: .4e}  analytic 18M^3H^2D={row['analytic_18M^3H^2D']: .4e}  "
              f"(single-h0 formula {row['formula_with_single_h0']: .4e})  ratio code/analytic="
              f"{row['EFT_kinetic_code'] / row['analytic_18M^3H^2D']: .3f}")
    print(f"  code verdict: {r[tag]['verdict']}")
print("EXPECTED: code/analytic = 1 at all a; model B rejected (scalar ghost, D < 0 for all a).")

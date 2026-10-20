"""Repro (track 6): the EFT_minkowski_limit check evaluates the tachyonic condition with EFTLambdadotdot = 0 for
full-mapping models that do not compute it analytically (quintessence 5e, K-mouflage, ADE, Horava, Horndeski, Extended
Galileon, OL alpha, designer f(R)/mc5e ...).

fortran/eftcamb/09_EFTCAMB_stability.f90:568-575 (EFTStabilityComputation):
    call input_model%compute_stability_factors(...)
    call input_model%compute_positivity_bounds(...)            <- uses EFTLambdadotdot, EFTGamma2PP(P), EFTGamma3PP(PP)
    if (a > pert_turn_on .and. EFT_mass_stability) then
        call input_model%compute_additional_derivs(...)         <- the only place these are filled for such models
EFTTestStability calls eft_cache%initialize() first, so EFTLambdadotdot is exactly 0 when the positivity / Minkowski
coefficients are built, whatever EFT_mass_stability is.

For quintessence (gamma_i = 0) the code's own formula reduces to
    ghostpos = D1 = -2 EFTc,   tachpos = -Q1/D1,   Q1 = (EFTLambdadotdot - adotoa EFTLambdadot) a^2,
and the Minkowski check fails where -ghostpos <= 0 or -tachpos < 0.  We rebuild EFTLambdadotdot from a centred finite
difference of the cached EFTLambdadot (EFTCAMB convention X = Lambda a^2: EFTLambdadotdot = d(EFTLambdadot)/dtau
- 2 adotoa EFTLambdadot) and evaluate the same condition.  The check interval is set to [1e-3, 1] (EFTCAMB_stability_time = 1e-3) because
for a < 1e-5 the frozen-field background makes EFTc ~ 1e-50 and the signs are roundoff.

Run: /Users/marco/Software/miniforge3/bin/python minkowski-tachpos-lambdadotdot-zero.py [tree]
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
PAR = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"

WORKER = r"""
import sys, json, ctypes
from ctypes import byref, c_double, c_int
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
from camb.baseconfig import camblib as lib
from camb.eftcamb import EFTCAMB_timestep_cache as C

out = {}
ini = sys.argv[2]
for mass in (False, True):
    cp = camb.read_ini(ini, no_validate=True)
    cp.EFTCAMB.EFT_minkowski_limit = True
    cp.EFTCAMB.EFT_mass_stability = mass
    try:
        camb.get_background(cp); verdict = "stable"
    except Exception as e:
        verdict = "unstable"
    cp = camb.read_ini(ini, no_validate=True)
    cp.EFTCAMB.EFT_minkowski_limit = True
    cp.EFTCAMB.EFT_mass_stability = mass
    cp.EFTCAMB.EFTCAMB_skip_stability = True
    P = camb.get_background(cp).Params
    test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
    buf = (ctypes.c_char * (1 << 18))()
    a = np.logspace(-3, 0, 3001)
    F = ("adotoa", "EFTc", "EFTLambdadot", "EFTLambdadotdot", "ghostpos", "tachpos", "EFTGamma2V", "EFTGamma3V")
    v = {f: np.empty_like(a) for f in F}
    ok = np.empty(len(a), bool)
    for i, x in enumerate(a):
        ok[i] = bool(test(byref(c_double(x)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf))
        raw = bytes(buf[:8 * 200])
        for f in F:
            v[f][i] = np.frombuffer(raw, np.float64, 1, getattr(C, f).offset)[0]
    H, Ld, D1 = v["adotoa"], v["EFTLambdadot"], v["ghostpos"]
    Ldd = H * np.gradient(Ld, np.log(a), edge_order=2) - 2 * H * Ld
    tach = -(a**2) * (Ldd - H * Ld) / D1
    s = slice(5, -5)
    fail_code = (-D1 <= 0) | (-v["tachpos"] < 0)
    fail_num = (-D1 <= 0) | (-tach < 0)
    i5 = int(np.argmin(abs(a - 0.5)))
    out[str(mass)] = {"verdict": verdict, "per_point_unstable_code": int((~ok)[s].sum()), "n": int(len(a[s])),
                      "max|gamma2|+|gamma3|": float(abs(v["EFTGamma2V"]).max() + abs(v["EFTGamma3V"]).max()),
                      "max|EFTLambdadotdot| used": float(abs(v["EFTLambdadotdot"]).max()),
                      "minkowski_fail_code": int(fail_code[s].sum()), "minkowski_fail_with_fd_Ldd": int(fail_num[s].sum()),
                      "a=0.5": {"EFTLambdadotdot_code": float(v["EFTLambdadotdot"][i5]), "EFTLambdadotdot_fd": float(Ldd[i5]),
                                "tachpos_code": float(v["tachpos"][i5]), "tachpos_fd": float(tach[i5])}}
print("RESULT " + json.dumps(out))
"""

with tempfile.TemporaryDirectory() as d:
    rows = {}
    for model in ("5_quint_exponential", "5_quint_inverse_power"):
        ini = Path(d) / f"{model}.ini"
        ini.write_text(f"DEFAULT({PAR}/{model}.ini)\nfeedback_level = 0\nEFTCAMB_write_background = F\n"
                       f"output_root = {d}/out\nEFTCAMB_stability_time = 1.e-3\n")
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        p = _run([sys.executable, "-c", WORKER, TREE, str(ini)], capture_output=True, text=True, cwd=d,
                           env=env)
        res = [json.loads(l[7:]) for l in p.stdout.splitlines() if l.startswith("RESULT ")]
        if not res:
            print(p.stdout[-2000:], p.stderr[-2000:])
            sys.exit(1)
        rows[model] = res[0]
print("tree:", TREE)
for model, r in rows.items():
    for mass, x in r.items():
        print(f"{model} EFT_mass_stability={mass}: code verdict (minkowski on) = {x['verdict']}; "
              f"max|gamma2|+|gamma3| = {x['max|gamma2|+|gamma3|']:.1e}; max|EFTLambdadotdot| used = "
              f"{x['max|EFTLambdadotdot| used']:.3e}")
        print(f"    Minkowski condition fails at {x['minkowski_fail_code']}/{x['n']} points with the code's coefficients, "
              f"at {x['minkowski_fail_with_fd_Ldd']}/{x['n']} with the finite-difference EFTLambdadotdot")
        y = x["a=0.5"]
        print(f"    a=0.5: EFTLambdadotdot code={y['EFTLambdadotdot_code']:.4e} fd={y['EFTLambdadotdot_fd']:.4e}; "
              f"tachpos code={y['tachpos_code']:.4e} fd={y['tachpos_fd']:.4e}")
print("EXPECTED: the Minkowski verdict is computed with the model's actual EFTLambdadotdot (here nonzero), i.e. the "
      "fail counts in the two columns agree; OBSERVED: EFTLambdadotdot = 0 is used and the verdict flips.")

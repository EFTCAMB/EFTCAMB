"""Repro (track 6): the EFT_mass_stability verdict is not reliable.

(a) Designer models (Pure EFT, RPH, ...) compute EFTcdotdot, EFTGamma1PP, EFTGamma2PP analytically, but
    EFTStabilityComputation (09_EFTCAMB_stability.f90:572-575) calls compute_additional_derivs before
    compute_mass_factors, and the generic designer implementation (06p2_abstract_EFTCAMB_designer.f90:72-79) overwrites
    them with finite differences of helpers that do not depend on the step (see repro additional-derivs-wrong-cache.py),
    i.e. EFTGamma1PP = EFTGamma2PP = 0 and EFTcdotdot = -2 adotoa EFTcdot.  The mass verdict of the test inis that
    enable EFT_mass_stability is computed from these values.
(b) With the check on, canonical quintessence (5_quint_exponential) is rejected: the code's mass eigenvalues reach
    |mu|/H^2 ~ 1e9 while the decoupling-limit pi mass of a canonical scalar, m^2 = -3 dH/dt - (3/2) H cdot/c, is O(H^2).

Part (a) compares the cache of the same model with EFT_mass_stability = F (model's own analytic values) and = T
(after compute_additional_derivs) at a = 0.5, against a centred finite difference of the first derivatives.
If the track-6 candidate-fix tree (fixaddl: helpers evaluated at temp_a in temp_cache) exists, the minimum mu/H^2 and
the verdict are also shown with correct derivatives.

Run: /Users/marco/Software/miniforge3/bin/python mass-stability-verdict-unreliable.py [tree]
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
FIX = f"{ROOT}/work/t6-stability/fixaddl"
PAR = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"

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
cp = camb.read_ini(ini, no_validate=True)
cp.EFTCAMB.EFT_mass_stability = True
try:
    camb.get_background(cp); out["verdict_mass_on"] = "stable"
except Exception:
    out["verdict_mass_on"] = "unstable"
test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
F = ("a", "adotoa", "Hdot", "EFTc", "EFTcdot", "EFTcdotdot", "EFTGamma1P", "EFTGamma1PP", "EFTGamma2P", "EFTGamma2PP",
     "EFT_mu1", "EFT_mu2")
def sample(P, grid):
    buf = (ctypes.c_char * (1 << 18))()
    v = {f: np.empty(len(grid)) for f in F}
    for i, a in enumerate(grid):
        test(byref(c_double(a)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf)
        raw = bytes(buf[:8 * 200])
        for f in F:
            v[f][i] = np.frombuffer(raw, np.float64, 1, getattr(C, f).offset)[0]
    return v
for mass in (False, True):
    cp = camb.read_ini(ini, no_validate=True)
    cp.EFTCAMB.EFTCAMB_skip_stability = True
    P = camb.get_background(cp).Params
    P.EFTCAMB.EFT_mass_stability = mass
    eps = 1e-4
    v = sample(P, np.array([0.5 * (1 - eps), 0.5, 0.5 * (1 + eps)]))
    H = v["adotoa"][1]
    fd = {"EFTGamma1PP": (v["EFTGamma1P"][2] - v["EFTGamma1P"][0]) / (eps),
          "EFTGamma2PP": (v["EFTGamma2P"][2] - v["EFTGamma2P"][0]) / (eps),
          "EFTcdotdot": 0.5 * H * (v["EFTcdot"][2] - v["EFTcdot"][0]) / eps - 2 * H * v["EFTcdot"][1]}
    out[f"mass_{mass}"] = {f: float(v[f][1]) for f in ("EFTcdotdot", "EFTGamma1PP", "EFTGamma2PP")}
    out["finite_difference"] = {k: float(x) for k, x in fd.items()}
    if mass:
        grid = np.logspace(-2, 0, 201)
        w = sample(P, grid)
        sel = grid > float(P.EFTCAMB.EFTCAMB_pert_turn_on)
        x = np.minimum(w["EFT_mu1"], w["EFT_mu2"])[sel] * (grid[sel] / w["adotoa"][sel]) ** 2
        mpi = (-3 * (w["Hdot"] / w["adotoa"] ** 2 - 1) - 1.5 * w["EFTcdot"] / (w["adotoa"] * w["EFTc"]))[sel]
        out["min_mu_over_H2"] = float(np.nanmin(x))
        out["max_abs_mu_over_H2"] = float(np.nanmax(np.abs(x)))
        out["decoupling_mpi2_over_H2_range"] = [float(np.nanmin(mpi)), float(np.nanmax(mpi))]
        out["rate"] = float(P.EFTCAMB.EFT_mass_stability_rate)
print("RESULT " + json.dumps(out))
"""


def run(tree, model):
    with tempfile.TemporaryDirectory() as d:
        ini = Path(d) / "m.ini"
        ini.write_text(f"DEFAULT({PAR}/{model}.ini)\nfeedback_level = 0\nEFTCAMB_write_background = F\n"
                       f"output_root = {d}/out\n")
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        p = _run([sys.executable, "-c", WORKER, tree, str(ini)], capture_output=True, text=True, cwd=d,
                           env=env, timeout=600)
    res = [json.loads(l[7:]) for l in p.stdout.splitlines() if l.startswith("RESULT ")]
    if not res:
        print(p.stdout[-1500:], p.stderr[-1500:])
        sys.exit(1)
    return res[0]


trees = [TREE] + ([FIX] if Path(FIX, "camb", "camblib.so").exists() else [])
print("trees:", trees)
for model in ("3_RPH_MB_const_1", "3_RPH_MK_const_1", "2_PEFT_Omega_lin_1", "5_quint_exponential"):
    for tree in trees:
        r = run(tree, model)
        tag = "A0/tree" if tree == TREE else "fixaddl"
        print(f"{model} [{tag}] verdict with EFT_mass_stability=T: {r['verdict_mass_on']}; min mu/H^2 after turn-on = "
              f"{r['min_mu_over_H2']:.4g} (threshold -rate^2 = {-r['rate']**2:.0f}); max |mu|/H^2 = "
              f"{r['max_abs_mu_over_H2']:.3g}" + (f"; decoupling-limit m_pi^2/H^2 of the canonical field in "
              f"[{r['decoupling_mpi2_over_H2_range'][0]:.3g}, {r['decoupling_mpi2_over_H2_range'][1]:.3g}]"
              if model.startswith("5_quint") else ""))
        if tree == TREE:
            for f in ("EFTcdotdot", "EFTGamma1PP", "EFTGamma2PP"):
                print(f"    a=0.5 {f:12s} model's own (mass=F) = {r['mass_False'][f]: .5e}   used by mass check "
                      f"(mass=T) = {r['mass_True'][f]: .5e}   finite difference = {r['finite_difference'][f]: .5e}")
print("EXPECTED: the mass check uses the model's own (finite-difference-consistent) derivatives, and a canonical "
      "quintessence field (|m_pi^2| ~ H^2) is not rejected with rate 10.")

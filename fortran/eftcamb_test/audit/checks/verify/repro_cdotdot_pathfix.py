"""Repro (track 6): the stability verdict depends on uninitialised memory.

EFTCAMB_Stability_Check (09_EFTCAMB_stability.f90:71) declares a local `type(TEFTCAMB_timestep_cache) :: eft_cache`
(no default initialisation).  Every EFTTestStability call does `call eft_cache%initialize()`, but
EFTCAMBTimestepCacheInit (03_abstract_EFTCAMB_cache.f90:402-586) resets neither EFTcdotdot nor EFTGamma1PP, while
EFTCAMBTimestepCacheIsNan checks both (lines 631, 637) and a NaN makes the model "unstable".  For models that never
assign them (all full-mapping models, designer f(R), Horava, ...; they are only filled by compute_additional_derivs,
i.e. with EFT_mass_stability = T and a > turn-on) the verdict is decided by whatever bytes were on the stack.

Part 1 (A0, deterministic): call the public EFTTestStability on the same model and scale factor with the cache
argument pre-filled with zeros vs pre-filled with NaN bit patterns (what -finit-real=nan -finit-derived, ifort -init=snan
or plain stack garbage can give).
Part 2 (optional): if the track-6 build with -finit-real=nan -finit-derived exists, run get_background with the stability
check on for a few test models.

Run: /Users/marco/Software/miniforge3/bin/python stability-cache-cdotdot-uninitialised.py [tree]
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
NANTREE = f"{ROOT}/work/t6-stability/nanini"
PAR = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"

WORKER = r"""
import sys, json, ctypes, struct
from ctypes import byref, c_double, c_int
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
from camb.baseconfig import camblib as lib
from camb.eftcamb import EFTCAMB_timestep_cache as C

mode = sys.argv[3]
cp = camb.read_ini(sys.argv[2], no_validate=True)
out = {}
if mode == "ctypes":
    cp.EFTCAMB.EFTCAMB_skip_stability = True
    P = camb.get_background(cp).Params
    test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
    for fill in ("zeros", "nan"):
        n = 1 << 15
        buf = (ctypes.c_double * n)(*([0.0] * n if fill == "zeros" else [float("nan")] * n))
        ok = bool(test(byref(c_double(0.5)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache),
                       buf))
        raw = bytes(buf)
        g = lambda f: float(np.frombuffer(raw, np.float64, 1, getattr(C, f).offset)[0])
        out[fill] = {"stable": ok, "EFTcdotdot": g("EFTcdotdot"), "EFTGamma1PP": g("EFTGamma1PP"),
                     "EFTcdot": g("EFTcdot"), "EFT_kinetic": g("EFT_kinetic")}
else:
    try:
        camb.get_background(cp); out["verdict"] = "stable"
    except Exception as e:
        out["verdict"] = "unstable"
print("RESULT " + json.dumps(out))
"""


def run(tree, model, mode, extra=""):
    with tempfile.TemporaryDirectory() as d:
        ini = Path(d) / "m.ini"
        ini.write_text(f"DEFAULT({PAR}/{model}.ini)\nfeedback_level = 3\nEFTCAMB_write_background = F\n"
                       f"output_root = {d}/out\nEFTCAMB_skip_stability = F\n{extra}")
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        p = _run([sys.executable, "-c", WORKER, tree, str(ini), mode], capture_output=True, text=True,
                           cwd=d, env=env, timeout=600)
    res = [json.loads(l[7:]) for l in p.stdout.splitlines() if l.startswith("RESULT ")]
    nan_lines = [l.strip() for l in p.stdout.splitlines() if " is NaN" in l]
    if not res:
        print(p.stdout[-1500:], p.stderr[-1500:])
        sys.exit(1)
    return res[0], nan_lines[:1]


print("Part 1, tree:", TREE)
for model in ("5_quint_exponential", "4_Des_F_R_LCDM_1"):
    r, msg = run(TREE, model, "ctypes")
    for fill, x in r.items():
        print(f"  {model}: EFTTestStability(a=0.5) with cache pre-filled with {fill:5s}: "
              f"{'stable' if x['stable'] else 'UNSTABLE'}  (EFTcdotdot after the call = {x['EFTcdotdot']}, "
              f"EFTGamma1PP = {x['EFTGamma1PP']})")
    print(f"     first NaN message: {msg}")
if Path(NANTREE, "camb", "camblib.so").exists():
    print("Part 2, tree built with -finit-real=nan -finit-derived:", NANTREE)
    for model in ("5_quint_exponential", "4_Des_F_R_LCDM_1", "5_Kmouflage_1", "2_PEFT_gamma1_const_1"):
        r0, _ = run(TREE, model, "verdict")
        r1, msg = run(NANTREE, model, "verdict")
        print(f"  {model:22s} release A0: {r0['verdict']:8s}  NaN-initialised build: {r1['verdict']:8s} {msg}")
print("EXPECTED: the verdict does not depend on the prior content of the cache (stable in every row).")

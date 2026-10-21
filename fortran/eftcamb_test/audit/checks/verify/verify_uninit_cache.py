"""Independent check of finding uninit-cache-fields-garbage (no rebuild, unmodified release library).

Check 1 (stab): call EFTTestStability (09_EFTCAMB_stability.f90) directly through ctypes for a model that CAMB already
    accepted (get_background succeeded), passing an eft_cache buffer whose memory is pre-filled with a chosen pattern
    (zeros / quiet NaN / NaN except EFTcdotdot+EFTGamma1PP). This emulates the uninitialised local
    `type(TEFTCAMB_timestep_cache) :: eft_cache` of EFTCAMB_Stability_Check holding a given bit pattern. After the call
    the buffer is read back and every field that is still NaN is listed: these are exactly the fields that neither
    initialize() nor the stability computation writes.
Check 2 (scribble): run with macOS MallocPreScribble=1 (every malloc'd block is filled with 0xAA bytes) and print
    get_eft_functions V_scf for Horndeski models with Horndeski_designer_V = F. If the Vscf table is never written the
    output is exactly the 0xAA pattern (-3.7206620809969885e-103) instead of 0.

Usage: python verify_uninit_cache.py [tree]     (tree label under audit/trees or absolute path; default A0)
"""

from __future__ import annotations

import json
import os
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
PARAMS = HERE.parents[2] / "parameters"
CHEAP = {"get_transfer": "F", "get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "100",
         "feedback_level": "0", "EFTCAMB_feedback_level": "0", "EFTCAMB_write_background": "F"}
SCRIBBLE = struct.unpack("<d", b"\xaa" * 8)[0]

STAB_MODELS = ["5_quint_exponential", "4_Des_F_R_LCDM_1", "5_Kmouflage_3", "5_ADE_1", "5_scg", "5_Horava_eta_1",
               "5_Kmimic_1", "5_hdsk_TG_1", "2_PEFT_Omega_const_1", "3_RPH_MK_const_1", "1_EFT_GR"]
SCRIB_MODELS = ["5_hdsk_TG_1", "5_hdsk_TG_2", "5_hdsk_TG_3"]


def wrap_ini(model: str, workdir: str, extra: dict | None = None) -> str:
    ini = Path(workdir) / f"{model}_wrap.ini"
    kv = dict(CHEAP, **(extra or {}))
    ini.write_text("\n".join(f"{k} = {v}" for k, v in kv.items()) + f"\nDEFAULT({PARAMS / (model + '.ini')})\n")
    return str(ini)


def child_stab(tree: str, model: str, workdir: str) -> None:
    sys.path.insert(0, tree)
    import ctypes

    import numpy as np

    import camb
    from camb.baseconfig import camblib
    from camb.eftcamb import EFTCAMB_timestep_cache as Cache

    cp = camb.read_ini(wrap_ini(model, workdir))
    res = camb.get_background(cp)  # runs the EFTCAMB init + full stability check: the model is accepted
    E = res.Params.EFTCAMB
    P = res.Params.EFTCAMB_parameter_cache
    test = getattr(camblib, "__eftcamb_stability_MOD_eftteststability")
    test.restype = ctypes.c_int
    cleanup = getattr(camblib, "__eftcamb_stability_MOD_eftstability_cleanup")
    names = [f for f, _ in Cache._fields_]
    out = {"model": model, "EFTflag": int(E.EFTflag), "astart": E.EFTCAMB_stability_time,
           "pert_turn_on": E.EFTCAMB_pert_turn_on, "mass_stability": bool(E.EFT_mass_stability), "tests": []}
    for a in (E.EFTCAMB_stability_time, 0.5):
        for fill in ("zeros", "nan", "nan_except_2"):
            buf = Cache()
            arr = np.ctypeslib.as_array(ctypes.cast(ctypes.byref(buf), ctypes.POINTER(ctypes.c_double)),
                                        shape=(ctypes.sizeof(Cache) // 8,))
            arr[:] = 0.0 if fill == "zeros" else np.nan
            if fill == "nan_except_2":
                buf.EFTcdotdot = 0.0
                buf.EFTGamma1PP = 0.0
            cleanup()
            ok = test(ctypes.byref(ctypes.c_double(a)), ctypes.byref(ctypes.c_double(10.0)), ctypes.byref(E),
                      ctypes.byref(P), ctypes.byref(buf))
            still_nan = [n for n in names if np.isnan(getattr(buf, n))]
            out["tests"].append({"a": a, "fill": fill, "stable": bool(ok), "still_nan": still_nan,
                                 "V_scf": getattr(buf, "V_scf")})
    print("JSON" + json.dumps(out), flush=True)


def child_scribble(tree: str, model: str, workdir: str) -> None:
    sys.path.insert(0, tree)
    import numpy as np

    import camb

    ini = wrap_ini(model, workdir)
    a = np.array([1e-3, 1e-2, 0.1, 0.5, 1.0])
    out = {"model": model, "runs": []}
    for i in range(3):
        cp = camb.read_ini(ini)
        r = camb.get_background(cp) if i != 1 else camb.get_results(cp)
        _, arr = cp.EFTCAMB.get_eft_functions(r, a)
        out["runs"].append({"V_scf": [float(x) for x in arr["V_scf"]], "phi_scf": [float(x) for x in arr["phi_scf"]],
                            "H0": float(r.Params.H0)})
    print("JSON" + json.dumps(out), flush=True)


def run_child(tree: str, mode: str, model: str, scribble: bool) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    if scribble:
        env["MallocPreScribble"] = "1"
    with tempfile.TemporaryDirectory() as d, harness.slot(1):
        p = subprocess.run(["nice", "-n", "10", sys.executable, __file__, "_child", mode, tree, model, d], env=env,
                           cwd=d, capture_output=True, text=True, timeout=900)
    for ln in p.stdout.splitlines():
        if ln.startswith("JSON"):
            return json.loads(ln[4:])
    return {"model": model, "error": f"rc={p.returncode}", "tail": (p.stdout + p.stderr)[-1500:]}


def main() -> None:
    arg = sys.argv[1] if len(sys.argv) > 1 else "A0"
    tree = arg if os.path.isabs(arg) else str(ROOT / "trees" / arg)
    which = sys.argv[2] if len(sys.argv) > 2 else "all"
    print(f"tree: {tree}")
    if which in ("all", "stab"):
        print("== Check 1: EFTTestStability with pre-filled eft_cache memory (unmodified library)")
        for m in STAB_MODELS:
            r = run_child(tree, "stab", m, scribble=False)
            if "error" in r:
                print(f"  {m}: child error {r['error']}\n{r['tail']}")
                continue
            print(f"  {m}: EFTflag={r['EFTflag']} astart={r['astart']:.3g} pert_turn_on={r['pert_turn_on']:.3g} "
                  f"mass_stability={r['mass_stability']}")
            for t in r["tests"]:
                print(f"     a={t['a']:<9.3g} fill={t['fill']:<13s} stable={t['stable']!s:<5s} "
                      f"fields still NaN after call={t['still_nan']}")
    if which in ("all", "scribble"):
        print(f"== Check 2: MallocPreScribble=1 (malloc fills with 0xAA = {SCRIBBLE!r}), get_eft_functions V_scf")
        for m in SCRIB_MODELS:
            r = run_child(tree, "scribble", m, scribble=True)
            if "error" in r:
                print(f"  {m}: child error {r['error']}\n{r['tail']}")
                continue
            for i, run in enumerate(r["runs"]):
                tag = ["get_background", "get_results", "get_background"][i]
                print(f"  {m} run {i + 1} ({tag}): V_scf={['%.6g' % x for x in run['V_scf']]}  "
                      f"all==0xAA pattern: {all(x == SCRIBBLE for x in run['V_scf'])}  "
                      f"phi_scf={['%.4g' % x for x in run['phi_scf']]}")
        print("  control without MallocPreScribble:")
        for m in SCRIB_MODELS[:1]:
            r = run_child(tree, "scribble", m, scribble=False)
            for i, run in enumerate(r.get("runs", [])):
                print(f"  {m} run {i + 1}: V_scf={['%.6g' % x for x in run['V_scf']]}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_child":
        mode, tree, model, d = sys.argv[2:6]
        (child_stab if mode == "stab" else child_scribble)(tree, model, d)
    else:
        main()

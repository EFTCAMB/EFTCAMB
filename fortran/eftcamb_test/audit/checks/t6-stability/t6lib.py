"""Track 6 (stability) probe library.

Runs, in a subprocess importing camb from a chosen frozen tree:
  1. camb.read_ini(DEFAULT(model.ini) + overrides) with EFTCAMB_skip_stability forced as requested,
  2. camb.get_background(cp)  (EFTCAMB background init + return-to-GR, the same path as get_results),
  3. direct ctypes calls of the Fortran stability routines on the State%CP objects:
       EFTTestStability(a, k_max, CP%EFTCAMB, CP%eft_par_cache, eft_cache)   -> verdict + full timestep cache
       EFTCAMB_Stability_Check(success, CP%EFTCAMB, CP%eft_par_cache, astart, aend, k_max)
     so the stability coefficients and the per-point verdict can be sampled on any grid.

Usage from other scripts:
    import t6lib
    res = t6lib.probe("A0", ini, overrides={...}, a_grid=[...], fields=[...], check=True)
    res = {"status", "bg_error", "full_check", "a", "verdict", <field>: array, "pert_turn_on", "log"}
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t6-stability"
WORK.mkdir(parents=True, exist_ok=True)

DEFAULT_FIELDS = [
    "a", "adotoa", "Hdot", "Hdotdot", "Hdotdotdot", "EFTOmegaV", "EFTOmegaP", "EFTOmegaPP", "EFTc", "EFTcdot",
    "EFTLambda", "EFTLambdadot", "EFTLambdadotdot", "EFTGamma1V", "EFTGamma2V", "EFTGamma2P", "EFTGamma2PP",
    "EFTGamma2PPP", "EFTGamma3V", "EFTGamma3P", "EFTGamma3PP", "EFTGamma3PPP", "EFTGamma3PPPP", "EFTGamma4V",
    "EFTGamma5V", "EFTGamma6V", "alphaB", "alphaK", "alphaT", "alphaM", "Meff2", "EFTpiA1", "EFTpiA2", "EFTpiB1",
    "EFTpiB2", "EFTpiC", "EFTpiD1", "EFTpiD2", "EFTAT", "EFTDT", "EFT_kinetic", "EFT_gradient", "EFT_mu1",
    "EFT_mu2", "posbound1", "posbound2", "ghostpos", "tachpos", "grhom_t", "gpresm_t", "grhob_t", "grhoc_t",
    "grhor_t", "grhog_t", "grhonu_tot", "gpinu_tot",
]


def _worker(cfg: dict) -> None:
    harness._exit_with_parent()
    sys.path.insert(0, cfg["tree"])
    import ctypes
    from ctypes import byref, c_double, c_int

    import camb
    from camb.eftcamb import EFTCAMB_timestep_cache

    assert Path(camb.__file__).resolve().is_relative_to(Path(cfg["tree"]).resolve()), camb.__file__
    out: dict = {"camb_file": camb.__file__}
    cp = camb.read_ini(cfg["ini"], no_validate=True)
    for path, value in cfg.get("python_set", {}).items():
        obj = cp
        parts = path.split(".")
        for p in parts[:-1]:
            obj = getattr(obj, p)
        setattr(obj, parts[-1], value)
    out["EFTflag"] = int(cp.EFTCAMB.EFTflag)
    try:
        results = camb.get_background(cp)
        out["bg_error"] = None
    except Exception as e:  # noqa: BLE001
        out["bg_error"] = repr(e)
        results = None
    if cfg.get("then_ini"):
        # process re-use check: run a second model in the same process after the first one
        cp2 = camb.read_ini(cfg["then_ini"], no_validate=True)
        try:
            r2 = camb.get_background(cp2)
            out["then_ok"] = True
            out["then_H0"] = float(r2.hubble_parameter(0.0))
        except Exception as e:  # noqa: BLE001
            out["then_ok"] = False
            out["then_error"] = repr(e)
    a_grid = np.asarray(cfg.get("a_grid", []), dtype=float)
    if results is not None and out["EFTflag"] != 0 and (len(a_grid) or cfg.get("check")):
        P = results.Params
        eft = P.EFTCAMB
        pc = P.EFTCAMB_parameter_cache
        out["pert_turn_on"] = float(eft.EFTCAMB_pert_turn_on)
        out["stability_time"] = float(eft.EFTCAMB_stability_time)
        for flag in ("EFT_ghost_math_stability", "EFT_mass_math_stability", "EFT_ghost_stability",
                     "EFT_gradient_stability", "EFT_mass_stability", "EFT_additional_priors",
                     "EFT_positivity_bounds", "EFT_minkowski_limit"):
            out[flag] = bool(getattr(eft, flag))
        out["EFT_mass_stability_rate"] = float(eft.EFT_mass_stability_rate)
        if "feedback" in cfg:
            eft.EFTCAMB_feedback_level = int(cfg["feedback"])
        from camb.baseconfig import camblib as lib
        test = getattr(lib, "__eftcamb_stability_MOD_eftteststability")
        test.restype = c_int
        cleanup = getattr(lib, "__eftcamb_stability_MOD_eftstability_cleanup")
        check = getattr(lib, "__eftcamb_stability_MOD_eftcamb_stability_check")
        buf = (ctypes.c_char * (1 << 18))()
        fields = cfg.get("fields") or DEFAULT_FIELDS
        offs = {f: getattr(EFTCAMB_timestep_cache, f).offset for f in fields}
        vals = {f: np.empty(len(a_grid)) for f in fields}
        verdict = np.empty(len(a_grid), dtype=bool)
        kmax = cfg.get("kmax", 10.0)
        cleanup()
        for i, a in enumerate(a_grid):
            if cfg.get("cleanup_each", True):
                cleanup()
            k = c_double(kmax)
            verdict[i] = bool(test(byref(c_double(a)), byref(k), byref(eft), byref(pc), buf))
            raw = bytes(buf[: max(offs.values()) + 8])
            for f, o in offs.items():
                vals[f][i] = np.frombuffer(raw, dtype=np.float64, count=1, offset=o)[0]
        out["a"] = a_grid.tolist()
        out["verdict"] = verdict.tolist()
        for f in fields:
            out[f] = vals[f].tolist()
        if cfg.get("check"):
            succ = c_int(0)
            k = c_double(kmax)
            astart = cfg.get("astart", out["stability_time"])
            aend = cfg.get("aend", 1.0)
            check(byref(succ), byref(eft), byref(pc), byref(c_double(astart)), byref(c_double(aend)), byref(k))
            out["full_check"] = bool(succ.value)
    print("T6RESULT " + json.dumps(out))


def probe(tree, ini, overrides=None, a_grid=(), fields=None, check=False, python_set=None, kmax=10.0,
          feedback=None, then_ini=None, timeout=900, name=None, astart=None, aend=None, threads=1,
          cleanup_each=True) -> dict:
    tree = harness.tree_path(tree)
    tmpdir = Path(tempfile.mkdtemp(prefix=(name or "probe") + "_", dir=WORK / "tmp" if (WORK / "tmp").exists()
                                   else WORK))
    ov = dict(overrides or {})
    lines = [f"DEFAULT({Path(ini).resolve()})", f"output_root = {tmpdir / 'out'}", "EFTCAMB_write_background = F"]
    lines += [f"{k} = {v}" for k, v in ov.items()]
    wrapper = tmpdir / "wrapper.ini"
    wrapper.write_text("\n".join(lines) + "\n")
    cfg = {"tree": str(tree), "ini": str(wrapper), "a_grid": list(map(float, a_grid)), "fields": fields,
           "check": check, "python_set": python_set or {}, "kmax": kmax, "cleanup_each": cleanup_each}
    if feedback is not None:
        cfg["feedback"] = feedback
    if then_ini:
        cfg["then_ini"] = str(Path(then_ini).resolve())
    if astart is not None:
        cfg["astart"] = astart
    if aend is not None:
        cfg["aend"] = aend
    cfg_path = tmpdir / "cfg.json"
    cfg_path.write_text(json.dumps(cfg))
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    env.pop("PYTHONPATH", None)
    try:
        with harness.slot(threads):
            proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), "_worker", str(cfg_path)],
                                  cwd=tmpdir, env=env, capture_output=True, text=True, timeout=timeout,
                                  preexec_fn=harness._nice)
        rc = proc.returncode
        stdout, stderr = proc.stdout, proc.stderr
    except subprocess.TimeoutExpired as e:
        rc, stdout, stderr = None, (e.stdout or b"").decode(errors="replace") if isinstance(e.stdout, bytes) else (
            e.stdout or ""), ""
    (tmpdir / "stdout.log").write_text(stdout + "\n--- stderr ---\n" + stderr)
    res: dict = {"status": "ok" if rc == 0 else ("timeout" if rc is None else f"rc={rc}"), "dir": str(tmpdir)}
    for line in stdout.splitlines():
        if line.startswith("T6RESULT "):
            res.update(json.loads(line[len("T6RESULT "):]))
    res["log"] = [l for l in stdout.splitlines() if not l.startswith("T6RESULT ")]
    res["stderr"] = stderr.splitlines()[-30:]
    for k in list(res):
        if isinstance(res[k], list) and res[k] and isinstance(res[k][0], (int, float, bool)) and k not in ("log",):
            res[k] = np.asarray(res[k])
    return res


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "_worker":
        _worker(json.loads(Path(sys.argv[2]).read_text()))

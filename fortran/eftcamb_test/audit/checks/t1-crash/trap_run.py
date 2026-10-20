"""Run harness models with gfortran floating-point traps ENABLED inside the Python process.

gfortran's -ffpe-trap only takes effect when the Fortran main program calls _gfortran_set_fpe at start-up. When CAMB
is used as a shared library from Python there is no Fortran main, so -ffpe-trap in the AD/ADD builds does nothing.
This runner calls _gfortran_set_fpe(13) (invalid|zero|overflow) and _gfortran_set_options (backtrace handlers) via
ctypes, but only while Fortran is running (camb.get_results and the EFTCAMB getters), so numpy code is unaffected.
On this machine (Apple arm64) a trapped FP exception is delivered as SIGILL.

Usage (same arguments as run_checked.py):
    python trap_run.py TREE OUTDIR [--models a,b | --pattern '5_*.ini' | --inis x.ini,y.ini] [--jobs 4]
    python trap_run.py --summary OUTDIR
Must use threads=1: the FP environment is per thread.
"""

from __future__ import annotations

import argparse
import ctypes
import json
import os
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

LIBGFORTRAN = "/Users/marco/Software/miniforge3/lib/libgfortran.5.dylib"
FPE_FLAGS = 1 | 4 | 8  # GFC_FPE_INVALID | GFC_FPE_ZERO | GFC_FPE_OVERFLOW


class _DlInfo(ctypes.Structure):
    _fields_ = [("dli_fname", ctypes.c_char_p), ("dli_fbase", ctypes.c_void_p), ("dli_sname", ctypes.c_char_p),
                ("dli_saddr", ctypes.c_void_p)]


def _print_image_bases(camb) -> None:
    """Print load addresses of camblib and libgfortran so raw backtrace addresses can be resolved with atos."""
    libc = ctypes.CDLL(None)
    from camb import baseconfig

    lib = baseconfig.camblib
    for handle, sym in ((lib, "__handles_MOD___copy_handles_C_cltransferdata"),
                        (ctypes.CDLL(LIBGFORTRAN), "_gfortran_set_fpe")):
        try:
            addr = ctypes.cast(getattr(handle, sym), ctypes.c_void_p).value
        except AttributeError:
            continue
        info = _DlInfo()
        if libc.dladdr(ctypes.c_void_p(addr), ctypes.byref(info)):
            print(f"IMAGE_BASE {info.dli_fname.decode()} {hex(info.dli_fbase)}", flush=True)


def _install_traps(tree: str) -> None:
    sys.path.insert(0, tree)
    import camb
    import camb.eftcamb as ce

    _print_image_bases(camb)
    lib = ctypes.CDLL(LIBGFORTRAN)
    opts = (ctypes.c_int * 7)(10308, 16383, 0, 1, 1, 1, 31)  # what gfortran -fcheck=all -fbacktrace main passes
    lib._gfortran_set_options(7, opts)

    def wrap(fn):
        def inner(*a, **k):
            lib._gfortran_set_fpe(FPE_FLAGS)
            try:
                return fn(*a, **k)
            finally:
                lib._gfortran_set_fpe(0)

        return inner

    camb.get_results = wrap(camb.get_results)
    ce.EFTCAMB.get_eft_functions = wrap(ce.EFTCAMB.get_eft_functions)
    ce.EFTCAMB.get_scale_evolution = wrap(ce.EFTCAMB.get_scale_evolution)


def run_model_trap(tree, ini, out, overrides=None, evol_k=None, timeout=7200, python_set=None) -> dict:
    out = Path(out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    ini = Path(ini)
    wrapper = harness._wrapper_ini(ini, out.parent, out.stem, overrides or {})
    cfg = {"tree": str(harness.tree_path(tree)), "ini": str(wrapper), "model_ini": str(ini.resolve()), "out": str(out),
           "evol_k": evol_k or [], "python_set": python_set or {}}
    env = dict(os.environ, OMP_NUM_THREADS="1", PYTHONHASHSEED="0")
    env.pop("PYTHONPATH", None)
    log = out.with_suffix(".log")
    t0 = time.time()
    try:
        with harness.slot(1), open(log, "w") as f:
            t0 = time.time()
            proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), "_worker", json.dumps(cfg)],
                                  cwd=out.parent, env=env, stdout=f, stderr=subprocess.STDOUT, timeout=timeout,
                                  preexec_fn=harness._nice)
        rc = proc.returncode
        status = "ok" if rc == 0 else (f"signal:{signal.Signals(-rc).name}" if rc < 0 else f"error:rc={rc}")
    except subprocess.TimeoutExpired:
        rc, status = None, "timeout"
    info = {"tree": str(tree), "model_ini": str(ini), "overrides": overrides or {}, "python_set": python_set or {},
            "threads": 1, "fpe_traps": True, "status": status, "returncode": rc, "elapsed": round(time.time() - t0, 2),
            "out": str(out), "log": str(log)}
    if status != "ok":
        info["log_tail"] = log.read_text(errors="replace").splitlines()[-40:]
    out.with_suffix(".json").write_text(json.dumps(info, indent=2))
    return info


def main() -> None:
    if len(sys.argv) > 2 and sys.argv[1] == "_worker":
        cfg = json.loads(sys.argv[2])
        _install_traps(cfg["tree"])
        harness._worker(cfg)
        return
    p = argparse.ArgumentParser()
    p.add_argument("tree", nargs="?")
    p.add_argument("outdir")
    p.add_argument("--pattern", default="*.ini")
    p.add_argument("--models", default="")
    p.add_argument("--inis", default="")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--timeout", type=float, default=7200)
    p.add_argument("--set", nargs="*", default=[])
    p.add_argument("--evol-k", default="")
    p.add_argument("--summary", action="store_true")
    a = p.parse_args()
    if a.summary:
        sys.argv = [sys.argv[0], "--summary", a.outdir]
        import run_checked

        run_checked.summary(Path(a.outdir))
        return
    if a.inis:
        inis = [Path(x) for x in a.inis.split(",")]
    else:
        inis = harness.model_inis(a.pattern)
        if a.models:
            want = set(a.models.split(","))
            inis = [i for i in inis if i.stem in want]
    over = dict(s.split("=", 1) for s in a.set)
    evk = [float(x) for x in a.evol_k.split(",") if x]
    outdir = Path(a.outdir)
    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(lambda i: run_model_trap(a.tree, i, outdir / f"{i.stem}.npz", over, evk, a.timeout), inis))
    for r in res:
        print(f"{Path(r['model_ini']).stem:40s} {r['status']:16s} {r['elapsed']:8.1f}s")


if __name__ == "__main__":
    main()

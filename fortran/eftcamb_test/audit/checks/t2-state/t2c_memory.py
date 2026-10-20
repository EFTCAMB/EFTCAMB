"""Track 2c: memory growth over repeated EFTCAMB calculations in one process.

Each mode runs N iterations of (build params, get_results) in a fresh subprocess (OMP_NUM_THREADS=1) and records the
resident set size every 10 iterations (psutil).  Growth per iteration is the slope of a linear fit of RSS against
iteration over the second half of the run (after allocator warm-up).

Modes:
  gr_readini      read_ini(GR) each iteration                          (baseline)
  gr_setparams    set_params(GR) each iteration                        (baseline)
  eft_readini     read_ini each iteration, alternating 3 EFT models
  eft_setparams   set_params each iteration, alternating 3 EFT models
  eft_reinit      one CAMBparams, cp.EFTCAMB.initialize_parameters each iteration, alternating 3 EFT models
  eft_samecp      one CAMBparams, initialize_parameters with the same model each iteration
  eft_background  set_params + get_background only, alternating models (isolates SetParams/EFT init)

Usage: python t2c_memory.py MODE [N]      (prints JSON with the RSS samples and the fitted slope)
       python t2c_memory.py all [N]       (runs every mode, at most 3 at a time)
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t2-state" / "memory"
TREE = os.environ.get("T2C_TREE", "A0")
MODELS = ["2_PEFT_Omega_const_1", "4_Des_F_R_LCDM_1", "5_quint_exponential"]
CHEAP_INI = {
    "l_max_scalar": os.environ.get("T2C_LMAX", "100"), "get_tensor_cls": "F", "get_transfer": "F", "do_lensing": "F",
    "feedback_level": "0",  # (EFTCAMB_write_background = F is added by harness._wrapper_ini)
}
CHEAP_PY = dict(H0=70.0, ombh2=0.0226, omch2=0.112, mnu=0.06, tau=0.09, As=2.1e-9, ns=0.96,
                lmax=int(os.environ.get("T2C_LMAX", "100")), WantTransfer=False, DoLensing=False, feedback_level=0)
MODES = ["gr_readini", "gr_setparams", "eft_readini", "eft_setparams", "eft_reinit", "eft_samecp", "eft_background",
         "gr_background", "eft_readini_bg", "eft_reinit_bg", "eft_evol", "gr_readini_bg", "eft_readini_only",
         "gr_readini_only", "eft_setparams_only", "eft_reinit_only"]


def eft_keys(model: str) -> dict[str, str]:
    from t2b_state import eft_keys as k

    return k(model)


def malloc_stats() -> tuple[int, int]:
    """(blocks_in_use, size_in_use) summed over all malloc zones (macOS). Unlike RSS this is not affected by
    paging / memory compression on a loaded machine, and it sees every gfortran ALLOCATE."""
    import ctypes

    class MS(ctypes.Structure):
        _fields_ = [("blocks_in_use", ctypes.c_uint), ("size_in_use", ctypes.c_size_t),
                    ("max_size_in_use", ctypes.c_size_t), ("size_allocated", ctypes.c_size_t)]

    libc = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
    libc.malloc_zone_statistics.argtypes = [ctypes.c_void_p, ctypes.POINTER(MS)]
    st = MS()
    libc.malloc_zone_statistics(None, ctypes.byref(st))
    return int(st.blocks_in_use), int(st.size_in_use)


def child(mode: str, n: int, out: Path) -> None:
    sys.path.insert(0, str(harness.tree_path(TREE)))
    sys.path.insert(0, str(HERE))
    import psutil

    import camb

    proc = psutil.Process()
    out.parent.mkdir(parents=True, exist_ok=True)
    wrappers = {}
    for m in ["1_EFT_GR", *MODELS]:
        wrappers[m] = str(harness._wrapper_ini(harness.PARAMS_DIR / f"{m}.ini", out.parent, f"{mode}_{m}", CHEAP_INI))
    keys = {m: eft_keys(m) for m in ["1_EFT_GR", *MODELS]}

    def setp(m):
        params = dict(CHEAP_PY)
        params.update(keys[m])
        for _ in range(20):
            try:
                return camb.set_params(**params)
            except camb.baseconfig.CAMBUnknownArgumentError as e:
                params.pop(str(e).split(":")[-1].strip())
        raise RuntimeError

    import numpy as np

    a_grid = np.logspace(-4, 0, 50)
    cp_reuse = setp(MODELS[0]) if "reinit" in mode or "samecp" in mode else None
    samples, t0 = [], time.time()
    errors = 0
    for it in range(n):
        m = MODELS[it % len(MODELS)]
        try:
            if mode == "gr_readini":
                cp = camb.read_ini(wrappers["1_EFT_GR"])
                camb.get_results(cp)
            elif mode == "gr_setparams":
                cp = setp("1_EFT_GR")
                camb.get_results(cp)
            elif mode == "eft_readini":
                cp = camb.read_ini(wrappers[m])
                camb.get_results(cp)
            elif mode == "eft_setparams":
                cp = setp(m)
                camb.get_results(cp)
            elif mode == "eft_reinit":
                cp_reuse.EFTCAMB.initialize_parameters(cp_reuse, keys[m], print_header=False)
                camb.get_results(cp_reuse)
            elif mode == "eft_samecp":
                cp_reuse.EFTCAMB.initialize_parameters(cp_reuse, keys[MODELS[0]], print_header=False)
                camb.get_results(cp_reuse)
            elif mode == "eft_background":
                cp = setp(m)
                r = camb.get_background(cp)
                cp.EFTCAMB.get_eft_functions(r, a_grid)
            elif mode == "gr_background":
                cp = setp("1_EFT_GR")
                camb.get_background(cp)
            elif mode == "eft_readini_bg":
                cp = camb.read_ini(wrappers[m])
                camb.get_background(cp)
            elif mode == "gr_readini_bg":
                cp = camb.read_ini(wrappers["1_EFT_GR"])
                camb.get_background(cp)
            elif mode == "eft_readini_only":
                cp = camb.read_ini(wrappers[m])
            elif mode == "gr_readini_only":
                cp = camb.read_ini(wrappers["1_EFT_GR"])
            elif mode == "eft_setparams_only":
                cp = setp(m)
            elif mode == "eft_reinit_only":
                cp_reuse.EFTCAMB.initialize_parameters(cp_reuse, keys[m], print_header=False)
            elif mode == "eft_reinit_bg":
                cp_reuse.EFTCAMB.initialize_parameters(cp_reuse, keys[m], print_header=False)
                camb.get_background(cp_reuse)
            elif mode == "eft_evol":
                cp = setp(m)
                r = camb.get_results(cp)
                cp.EFTCAMB.get_scale_evolution(r, [0.01, 0.1], a_grid)
        except Exception as e:  # noqa: BLE001
            errors += 1
            if errors < 5:
                print(f"iteration {it}: {type(e).__name__}: {e}", flush=True)
        if it % int(os.environ.get('T2C_EVERY', '10')) == 0 or it == n - 1:
            import gc

            gc.collect()
            blocks, inuse = malloc_stats()
            samples.append((it, proc.memory_info().rss, time.time() - t0, inuse, blocks))
            print(f"{mode} it={it} rss={samples[-1][1] / 2**20:.1f}MiB malloc_in_use={inuse / 2**20:.3f}MiB "
                  f"blocks={blocks} t={samples[-1][2]:.0f}s", flush=True)
            out.write_text(json.dumps({"mode": mode, "n": n, "samples": samples, "errors": errors, "partial": True}))
    res = {"mode": mode, "n": n, "samples": samples, "errors": errors}
    if os.environ.get("T2C_LEAKS"):
        # macOS leaks: unreferenced malloc blocks in this (live) process; with MallocStackLogging=lite in the
        # environment the report carries the allocation backtrace of every leaked block.
        import gc
        import re

        gc.collect()
        p = subprocess.run(["leaks", str(os.getpid())], capture_output=True, text=True, timeout=400)
        out.with_suffix(".leaks.txt").write_text(p.stdout)
        m = re.search(r"(\d+) leaks for (\d+) total leaked bytes", p.stdout)
        res["leaks"] = [int(m.group(1)), int(m.group(2))] if m else None
        print(f"{mode} LEAKS n={n}: {res['leaks']}", flush=True)
    out.write_text(json.dumps(res))


def slope(samples: list, col: int = 1) -> tuple[float, float]:
    import numpy as np

    s = np.array(samples, dtype=float)
    half = s[s[:, 0] >= s[-1, 0] / 2]
    p = np.polyfit(half[:, 0], half[:, col], 1)
    return float(p[0]), float(s[-1, col] - s[0, col])


def main(modes: list[str], n: int) -> None:
    from concurrent.futures import ThreadPoolExecutor

    WORK.mkdir(parents=True, exist_ok=True)

    def one(mode):
        out = WORK / TREE / f"{mode}_{n}.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        nthr = int(env.get("T2C_THREADS", "1"))
        env["OMP_NUM_THREADS"] = str(nthr)
        with harness.slot(nthr), open(out.with_suffix(".log"), "w") as f:
            try:
                rc = subprocess.run([sys.executable, __file__, "_child", mode, str(n), str(out)], env=env, stdout=f,
                                    stderr=subprocess.STDOUT, cwd=out.parent, preexec_fn=harness._nice,
                                    timeout=float(os.environ.get("T2C_TIMEOUT", "520"))).returncode
            except subprocess.TimeoutExpired:
                rc = "timeout"
        return mode, rc, out

    with ThreadPoolExecutor(int(os.environ.get("T2C_JOBS", "3"))) as ex:
        res = list(ex.map(one, modes))
    for mode, rc, out in res:
        if not out.exists():
            print(f"{mode:16s} rc={rc} (see {out.with_suffix('.log')})")
            continue
        d = json.loads(out.read_text())
        if rc != 0 or d.get("partial"):
            print(f"{mode:16s} rc={rc} PARTIAL ({len(d['samples'])} samples)")
        if len(d["samples"]) < 3:
            continue
        sl, tot = slope(d["samples"])
        msl, mtot = slope(d["samples"], 3) if len(d["samples"][0]) > 3 else (float("nan"), float("nan"))
        bsl, btot = slope(d["samples"], 4) if len(d["samples"][0]) > 4 else (float("nan"), float("nan"))
        print(f"{mode:16s} n={d['n']} errors={d['errors']} rss0={d['samples'][0][1] / 2**20:.1f}MiB "
              f"rss_end={d['samples'][-1][1] / 2**20:.1f}MiB rss_growth={tot / 2**10:.0f}KiB "
              f"rss_slope(2nd half)={sl / 2**10:.2f}KiB/iter | malloc_in_use growth={mtot / 2**10:.1f}KiB "
              f"slope(2nd half)={msl / 2**10:.3f}KiB/iter blocks growth={btot:.0f} slope={bsl:.2f}/iter "
              f"time={d['samples'][-1][2]:.0f}s")


if __name__ == "__main__":
    if sys.argv[1] == "_child":
        harness._exit_with_parent()
        child(sys.argv[2], int(sys.argv[3]), Path(sys.argv[4]))
    else:
        nn = int(sys.argv[2]) if len(sys.argv) > 2 else 300
        main(MODES if sys.argv[1] == "all" else sys.argv[1].split(","), nn)

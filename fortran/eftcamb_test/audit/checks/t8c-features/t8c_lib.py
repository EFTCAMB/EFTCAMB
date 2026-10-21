"""Track 8c (CAMB v2 features with EFT on): subprocess runner for task-based CAMB worker scripts.

run_task(tree, ini, task, args, out, overrides) runs t8c_worker.py in a fresh python process that imports camb only
from the selected frozen tree, reads DEFAULT(ini)+overrides with camb.read_ini, applies the task (Python-side setup
and outputs) and writes out.npz (+ out.json status).
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

PARAMS = harness.PARAMS_DIR
INIS = HERE.parents[1] / "inis"
WORK = harness.AUDIT_ROOT / "work" / "t8c-features"
PY = sys.executable


def model_ini(name: str) -> Path:
    p = PARAMS / f"{name}.ini"
    if not p.exists():
        p = INIS / f"{name}.ini"
    return p


def run_task(tree, ini, task: str, out, args: dict | None = None, overrides: dict | None = None, threads: int = 1,
             timeout: float = 1800) -> dict:
    out = Path(out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    wrapper = harness._wrapper_ini(Path(ini), out.parent, out.stem, overrides or {})
    cfg = {"tree": str(harness.tree_path(tree)), "ini": str(wrapper), "task": task, "args": args or {},
           "out": str(out)}
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), PYTHONHASHSEED="0")
    env.pop("PYTHONPATH", None)
    log = out.with_suffix(".log")
    t0 = time.time()
    try:
        with harness.slot(threads), open(log, "w") as f:
            proc = subprocess.run([PY, str(HERE / "t8c_worker.py"), json.dumps(cfg)], cwd=out.parent, env=env,
                                  stdout=f, stderr=subprocess.STDOUT, timeout=timeout,
                                  preexec_fn=lambda: os.nice(10))
        rc = proc.returncode
        status = "ok" if rc == 0 else (f"signal:{signal.Signals(-rc).name}" if rc < 0 else f"error:rc={rc}")
    except subprocess.TimeoutExpired:
        rc, status = None, "timeout"
    info = {"tree": str(tree), "ini": str(ini), "task": task, "args": args or {}, "overrides": overrides or {},
            "status": status, "rc": rc, "elapsed": round(time.time() - t0, 1), "out": str(out), "log": str(log)}
    if status != "ok":
        info["log_tail"] = log.read_text(errors="replace").splitlines()[-30:]
    out.with_suffix(".json").write_text(json.dumps(info, indent=2))
    return info


def run_many(jobs: list[dict], nproc: int = 4) -> list[dict]:
    """jobs: list of kwargs for run_task."""
    with ThreadPoolExecutor(nproc) as ex:
        return list(ex.map(lambda kw: run_task(**kw), jobs))


def load(path) -> dict:
    return harness.load(path)


def maxrel(a, b, floor_frac: float = 0.0) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    den = np.abs(b)
    if floor_frac:
        den = np.maximum(den, floor_frac * np.max(np.abs(b)))
    m = den > 0
    return float(np.max(np.abs(a[m] - b[m]) / den[m])) if m.any() else 0.0


def cl_diffs(a, b, lmin=2, lmax=None) -> dict:
    """a, b: (L, 4) TT EE BB TE arrays. TE normalised by sqrt(TT EE)."""
    n = min(len(a), len(b)) if lmax is None else min(len(a), len(b), lmax + 1)
    a, b = np.asarray(a)[lmin:n], np.asarray(b)[lmin:n]
    res = {}
    for i, lab in enumerate(("TT", "EE", "BB")):
        if np.any(b[:, i] != 0):
            res[lab] = maxrel(a[:, i], b[:, i])
    if a.shape[1] > 3:
        norm = np.sqrt(np.abs(b[:, 0] * b[:, 1]))
        m = norm > 0
        res["TE"] = float(np.max(np.abs(a[m, 3] - b[m, 3]) / norm[m]))
    return res


def pk_diff(a: dict, b: dict, kmin: float = 1e-4, kmax: float | None = None, key: str = "pk",
            kkey: str = "pk_k") -> float:
    """Max fractional P(k) difference over common k (log-interpolated if grids differ), all redshifts."""
    ka, kb = a[kkey], b[kkey]
    if ka.shape == kb.shape and np.allclose(ka, kb, rtol=1e-12):
        sel = (ka >= kmin) & (ka <= (kmax or ka[-1]))
        return maxrel(a[key][:, sel], b[key][:, sel])
    lo, hi = max(ka[0], kb[0], kmin), min(ka[-1], kb[-1], kmax or np.inf)
    k = np.logspace(np.log10(lo), np.log10(hi), 500)
    out = 0.0
    for i in range(min(a[key].shape[0], b[key].shape[0])):
        pa = np.exp(np.interp(np.log(k), np.log(ka), np.log(a[key][i])))
        pb = np.exp(np.interp(np.log(k), np.log(kb), np.log(b[key][i])))
        out = max(out, maxrel(pa, pb))
    return out

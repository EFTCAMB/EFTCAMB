"""Independent verification of finding stability-sampler-stale-y (third stability pass never updates y).

Different model, window and methods from the finding's repro (which copies the loop formulas into numpy):

  A. "dump": the a values the sampler really visits, read from EFTCAMB's own DebugEFTCAMB cache dump
     (<outroot>stability_cache_FRW.dat, one row per EFTTestStability call) written by the ADD tree (-DDEBUG).
     Model C1 is stable for the A0 logic, so all three passes run to completion (3000 rows).
  B. "scan": on A0 (stability skipped), the scalar kinetic term is recomputed by hand from the EFT functions returned
     by the python bridge get_eft_functions (EFTCAMB_GetEFTFunctions, a code path that does not go through
     EFTTestStability) at every a the sampler visited (from A) and on a dense grid around the window.
  C. "consequence": full pipeline (harness.run_model, CMB+P(k)) for C1 (ghost condition) and C2 (only the
     A1-sign-change mathematical condition) on A0 and on the t6 fixed tree (FIXY = A0 + y update in loop 3).

Model: Pure EFT, LCDM background, Omega = 0, gamma1(a) = G0 + G2 (a - AW)^2 (Taylor series, no factorials),
ghost for |a - AW| < HW.  AW sits between two points of the linear pass and far from the log-near-astart pass.

Usage: python staley_indep.py [dump|scan|consequence|all]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

AUDIT = harness.AUDIT_ROOT
REPO = harness.REPO
BASE = REPO / "fortran/eftcamb_test/parameters/base_params.ini"
OUT = AUDIT / "work/verify/staley"
TREES = {"A0": AUDIT / "trees/A0", "ADD": AUDIT / "trees/ADD", "FIXY": AUDIT / "work/t6-stability/fixy"}

AW = 0.9974975  # midway between linear-pass points 996/999 and 997/999
HW = 3.0e-4
G0 = -1.0e-3
G2 = -G0 / HW**2


def ini_text(case: str, extra: str = "") -> str:
    if case == "C1":  # physical ghost condition only
        flags = "EFT_ghost_math_stability = F\nEFT_ghost_stability = T\nEFT_gradient_stability = T\n"
    else:  # C2: A1 sign change (mathematical) only
        flags = "EFT_ghost_math_stability = T\nEFT_ghost_stability = F\nEFT_gradient_stability = F\n"
    return f"""output_root = dbg
EFTflag = 1
PureEFTmodel = 1
PureEFTHorndeski = F
EFTwDE = 0
PureEFTmodelOmega = 0
PureEFTmodelGamma1 = 5
PureEFTmodelGamma2 = 0
PureEFTmodelGamma3 = 0
PureEFTmodelGamma4 = 0
PureEFTmodelGamma5 = 0
PureEFTmodelGamma6 = 0
EFTGamma1_Taylor_order = 2
EFTGamma1a0 = {AW!r}
EFTGamma10 = {G0!r}
EFTGamma11 = 0.
EFTGamma12 = {G2!r}
{flags}EFT_mass_math_stability = F
EFT_mass_stability = F
EFT_additional_priors = F
EFTCAMB_stability_time = 1.e-10
feedback_level = 0
{extra}DEFAULT({BASE})
"""


def write_ini(case: str, name: str, extra: str = "") -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / name
    p.write_text(ini_text(case, extra))
    return p


def run_worker(tree: Path, code: str, args: list[str], cwd: Path, timeout: float = 600) -> subprocess.CompletedProcess:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    cwd.mkdir(parents=True, exist_ok=True)
    with harness.slot(1):
        return subprocess.run([sys.executable, "-c", code, str(tree), *args], cwd=cwd, env=env, capture_output=True,
                              text=True, timeout=timeout, preexec_fn=lambda: os.nice(10))


DUMP_WORKER = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
print("EFTflag", cp.EFTCAMB.EFTflag, flush=True)
try:
    camb.get_background(cp)
    print("VERDICT stable", flush=True)
except Exception as e:
    print("VERDICT exception", str(e).splitlines()[-1], flush=True)
"""


def part_dump() -> dict:
    ini = write_ini("C1", "C1_dump.ini")
    d = OUT / "dump_ADD"
    for f in d.glob("dbg_*"):
        f.unlink()
    p = run_worker(TREES["ADD"], DUMP_WORKER, [str(ini)], d)
    frw = d / "dbg_stability_cache_FRW.dat"
    res = {"rc": p.returncode, "stdout_tail": p.stdout[-600:], "stderr_tail": p.stderr[-600:],
           "dump_exists": frw.exists()}
    if not frw.exists():
        return res
    a = np.loadtxt(frw, usecols=0, ndmin=1)
    np.save(OUT / "visited_a_ADD.npy", a)
    res["n_rows"] = int(a.size)
    inwin = (a > AW - HW) & (a < AW + HW)
    passes = {}
    for i in range(3):
        s = a[1000 * i:1000 * (i + 1)]
        w = (s > AW - HW) & (s < AW + HW)
        passes[f"pass{i + 1}"] = {"n": int(s.size), "unique": int(np.unique(s).size), "min": float(s.min()),
                                  "max": float(s.max()), "n_in_window": int(w.sum()),
                                  "n_above_0.99": int((s > 0.99).sum())}
    res["passes"] = passes
    res["n_in_window_total"] = int(inwin.sum())
    # what the y-updated third pass would visit, for comparison
    i = np.arange(1000)
    y = -10 + i * 10 / 999
    intended = 1 + (1e-10 - 1) * 10**y
    res["intended_pass3_n_in_window"] = int(((intended > AW - HW) & (intended < AW + HW)).sum())
    return res


SCAN_WORKER = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
cp.EFTCAMB.EFTCAMB_skip_stability = True
res = camb.get_background(cp)
a = np.load(sys.argv[3])
fields, c = cp.EFTCAMB.get_eft_functions(res, a)
h0 = cp.H0 / 299792.458
g = lambda n: np.asarray(c[n], dtype=float)
F = 1 + g("EFTOmegaV") - g("EFTGamma4V")
kin = 9 * F * (4 * g("EFTc") * F + 3 * g("adotoa")**2 * g("EFTOmegaP")**2 * a**2
               + a**2 * h0 * (h0 * (3 * g("EFTGamma2V")**2 + 8 * g("EFTGamma1V") * F)
                              + 6 * g("adotoa") * g("EFTGamma2V") * g("EFTOmegaP")))
np.savez(sys.argv[4], a=a, kin=kin, gamma1=g("EFTGamma1V"), c=g("EFTc"), omega=g("EFTOmegaV"))
print("DONE", flush=True)
"""


def part_scan() -> dict:
    visited = np.load(OUT / "visited_a_ADD.npy")
    dense = np.linspace(0.9960, 0.9990, 30001)
    arr = np.concatenate([visited, dense])
    order = np.argsort(arr, kind="stable")
    np.save(OUT / "scan_a.npy", arr[order])
    ini = write_ini("C1", "C1_scan.ini")
    p = run_worker(TREES["A0"], SCAN_WORKER, [str(ini), str(OUT / "scan_a.npy"), str(OUT / "scan_out.npz")],
                   OUT / "scan_A0")
    if "DONE" not in p.stdout:
        return {"error": p.stdout[-800:] + p.stderr[-800:]}
    z = np.load(OUT / "scan_out.npz")
    kin_sorted = z["kin"]
    kin = np.empty_like(kin_sorted)
    kin[order] = kin_sorted
    g1 = np.empty_like(kin_sorted)
    g1[order] = z["gamma1"]
    kv, kd = kin[:visited.size], kin[visited.size:]
    neg = dense[kd < 0]
    res = {
        "visited_points_with_negative_kinetic": {f"pass{i + 1}": int((kv[1000 * i:1000 * (i + 1)] < 0).sum())
                                                 for i in range(3)},
        "visited_min_kinetic_above_0.99": float(kv[visited > 0.99].min()),
        "dense_negative_kinetic_points": int(neg.size),
        "dense_negative_range": [float(neg.min()), float(neg.max())] if neg.size else None,
        "kinetic_at_AW": float(kd[np.argmin(abs(dense - AW))]),
        "gamma1_at_AW": float(g1[visited.size + np.argmin(abs(dense - AW))]),
        "analytic_window": [AW - HW, AW + HW],
        "max_abs_EFTc": float(np.abs(z["c"]).max()),
        "max_abs_Omega": float(np.abs(z["omega"]).max()),
        "all_visited_kinetic_positive": bool((kv > 0).all()),
    }
    return res


def part_consequence() -> dict:
    res = {}
    for case in ("C1", "C2"):
        ini = write_ini(case, f"{case}_full.ini")
        for t in ("A0", "FIXY"):
            st = harness.run_model(TREES[t], ini, OUT / "full" / f"{case}_{t}.npz", threads=1, timeout=900)
            meta = {}
            npz = OUT / "full" / f"{case}_{t}.npz"
            if npz.exists():
                z = harness.load(npz)
                m = z.get("meta")
                meta = m if isinstance(m, dict) else (json.loads(str(m)) if m is not None else {})
            sig8 = None
            if npz.exists() and st.get("status") == "ok":
                sig8 = float(np.atleast_1d(harness.load(npz)["sigma8"])[-1])
            res[f"{case}|{t}"] = {"status": st.get("status"), "errors": meta.get("errors"), "sigma8_z0": sig8,
                                  "sigma8": meta.get("derived", {}).get("sigma8") if isinstance(meta.get("derived"),
                                                                                               dict) else None,
                                  "log_tail": (st.get("log_tail") or "")[-300:]}
    return res


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    out = {}
    if what in ("dump", "all"):
        out["dump"] = part_dump()
    if what in ("scan", "all"):
        out["scan"] = part_scan()
    if what in ("consequence", "all"):
        out["consequence"] = part_consequence()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"staley_indep_{what}.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str))

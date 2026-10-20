"""Track 2a: thread determinism. Run each model with OMP_NUM_THREADS=1 and =8 on a tree and compare outputs.

Usage:
    python t2a_threads.py run [--tree A0] [--models m1,m2] [--tag t1] [--threads 1] [--jobs 4]
    python t2a_threads.py compare [--tree A0] [--a t1] [--b t8] [--models ...]

Outputs go to $WORK/threads/<tree>/<model>_<tag>.npz.  The compare step prints, for each model, the worst relative
difference of Cls / P(k) / sigma8 / background, the worst EFT-function difference (max|diff|/max|field|), the worst
get_scale_evolution difference, and whether everything is bit identical.
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t2-state"
EVOL_K = [1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.3, 1.0]

MODELS = [
    "1_EFT_GR",
    "2_PEFT_Omega_const_1", "2_PEFT_Omega_exp_1", "2_PEFT_Omega_lin_1", "2_PEFT_Omega_power_1",
    "2_PEFT_gamma1_const_1", "2_PEFT_gamma3_exp_1", "2_PEFT_gamma3_lin_1", "2_PEFT_gamma4_power_1",
    "2_PEFT_gamma5_exp_1", "2_PEFT_gamma6_lin_1", "2_PEFT_wDE_wCDM_1", "2_PEFT_wDE_LCDM_1",
    "3_RPH_KT_OmegaDE_1", "3_RPH_KT_power_1", "3_RPH_MB_OmegaDE_1", "3_RPH_MK_const_1", "3_RPH_MK_power_1",
    "4_Des_F_R_LCDM_1", "4_Des_F_R_wCDM_2", "4_Des_mc5e_CPL_1", "4_Des_mc5e_taylor_1",
    "4_Des_mc5e_turning_point_1", "4_Des_mc5e_wCDM_1",
    "5_ADE_1", "5_Horava_eta_1", "5_Horava_eta_4", "5_Horava_SolSyst_1", "5_Kmimic_1", "5_Kmouflage_1",
    "5_Kmouflage_3", "5_hdsk_TG_1", "5_hdsk_TG_3", "5_quint_axion", "5_quint_cubic", "5_quint_exponential", "5_scg",
]


def ini_path(model: str) -> Path:
    p = harness.PARAMS_DIR / f"{model}.ini"
    if p.exists():
        return p
    p = HERE.parents[1] / "inis" / f"{model}.ini"
    if p.exists():
        return p
    raise FileNotFoundError(model)


def outdir(tree: str) -> Path:
    d = WORK / "threads" / Path(tree).name
    d.mkdir(parents=True, exist_ok=True)
    return d


def run(tree: str, models: list[str], tag: str, threads: int, jobs: int, evol: bool) -> None:
    d = outdir(tree)

    def one(m: str) -> dict:
        out = d / f"{m}_{tag}.npz"
        if out.exists():
            return {"model": m, "status": "cached"}
        info = harness.run_model(tree, ini_path(m), out, threads=threads, evol_k=EVOL_K if evol else None,
                                 timeout=3000)
        print(f"{m:32s} {tag:6s} {info['status']:14s} {info['elapsed']:8.1f}s", flush=True)
        return info

    with ThreadPoolExecutor(jobs) as ex:
        list(ex.map(one, models))


def rel_field(x: np.ndarray, y: np.ndarray) -> float:
    if x.shape != y.shape:
        return float("inf")
    fin = np.isfinite(x) & np.isfinite(y)
    if not np.array_equal(np.isfinite(x), np.isfinite(y)):
        return float("inf")
    scale = np.max(np.abs(y[fin])) if fin.any() else 0.0
    if scale == 0:
        return 0.0 if np.array_equal(x[fin], y[fin]) else float("inf")
    return float(np.max(np.abs(x[fin] - y[fin])) / scale)


def compare_runs(pa: Path, pb: Path) -> dict:
    A, B = harness.load(pa), harness.load(pb)
    c = harness.compare(A, B)
    res = {
        "worst_cl": c["worst_cl"],
        "worst_pk": c["worst_pk"],
        "sigma8": c.get("sigma8", 0.0),
        "bg_H": c.get("bg_H", 0.0),
        "bg_tau": c.get("bg_tau", 0.0),
        "bit_identical_obs": c["bit_identical"],
    }
    eft_keys = sorted(k for k in A if k.startswith("eft_") and k != "eft_a_grid")
    evol_keys = sorted(k for k in A if k.startswith("evol_") and k not in ("evol_k", "evol_a"))
    worst_eft, worst_eft_key, eft_exact = 0.0, "", True
    for k in eft_keys:
        if k not in B:
            continue
        x, y = A[k], B[k]
        eft_exact &= bool(np.array_equal(x, y, equal_nan=True))
        r = rel_field(x, y)
        if r > worst_eft:
            worst_eft, worst_eft_key = r, k
    res.update(worst_eft=worst_eft, worst_eft_key=worst_eft_key, eft_exact=eft_exact)
    worst_ev, worst_ev_key, ev_exact, ev_per_k = 0.0, "", True, []
    if evol_keys:
        # note: the timestep-cache fields "k" and "a" overwrite the harness keys evol_k / evol_a (shape (nk, na))
        nk = A[evol_keys[0]].shape[0]
        ev_per_k = [0.0] * nk
        for k in evol_keys:
            if k not in B:
                continue
            x, y = A[k], B[k]
            ev_exact &= bool(np.array_equal(x, y, equal_nan=True))
            r = rel_field(x, y)
            if r > worst_ev:
                worst_ev, worst_ev_key = r, k
            # evol arrays are [time, k] in ctypes order (times inner): shape (nk, na)
            for ik in range(nk):
                ev_per_k[ik] = max(ev_per_k[ik], rel_field(x[ik], y[ik]))
    res.update(worst_evol=worst_ev, worst_evol_key=worst_ev_key, evol_exact=ev_exact, evol_per_k=ev_per_k)
    res["all_exact"] = res["bit_identical_obs"] and eft_exact and ev_exact
    return res


def compare(tree: str, models: list[str], a: str, b: str, out_json: Path | None) -> dict:
    d = outdir(tree)
    table = {}
    for m in models:
        pa, pb = d / f"{m}_{a}.npz", d / f"{m}_{b}.npz"
        if not (pa.exists() and pb.exists()):
            print(f"{m:32s} missing")
            continue
        r = compare_runs(pa, pb)
        table[m] = r
        print(f"{m:32s} cl={r['worst_cl']:.2e} pk={r['worst_pk']:.2e} s8={r['sigma8']:.1e} H={r['bg_H']:.1e} "
              f"eft={r['worst_eft']:.1e} evol={r['worst_evol']:.1e}({r['worst_evol_key']}) exact={r['all_exact']}")
    if out_json:
        out_json.write_text(json.dumps(table, indent=1))
    return table


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=["run", "compare"])
    p.add_argument("--tree", default="A0")
    p.add_argument("--models", default="")
    p.add_argument("--tag", default="t1")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--no-evol", action="store_true")
    p.add_argument("--a", default="t1")
    p.add_argument("--b", default="t8")
    p.add_argument("--json", default="")
    args = p.parse_args()
    models = args.models.split(",") if args.models else MODELS
    if args.cmd == "run":
        run(args.tree, models, args.tag, args.threads, args.jobs, not args.no_evol)
    else:
        compare(args.tree, models, args.a, args.b, Path(args.json) if args.json else None)


if __name__ == "__main__":
    main()

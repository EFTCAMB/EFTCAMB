"""Track 6: for every model ini, sample EFTTestStability on (a) the exact grid the code's sampler uses and (b) a much
finer grid, and compare with the code's overall verdict (EFTCAMB_Stability_Check).  Reconstructs which condition
fails at each point from the cached stability coefficients.

python scan_all.py [--tree A0] [--pattern '*.ini'] [--jobs 4] [--extra-inis DIR] [--tag NAME] [--set k=v ...]
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import t6lib  # noqa: E402
import harness  # noqa: E402

AST, AEND, N = 1e-10, 1.0, 1000


def sampler_grid(astart=AST, aend=AEND, fixed=False):
    i = np.arange(N)
    lin = astart + i * (aend - astart) / (N - 1)
    y = -10.0 + i * 10.0 / (N - 1)
    log1 = astart + (aend - astart) * 10.0**y
    if fixed:
        log2 = aend + (astart - aend) * 10.0**y
    else:
        log2 = np.full(N, aend + (astart - aend) * 10.0**y[-1])  # stale y = 0 -> a = astart, 1000 times
    return lin, log1, log2


def fine_grid(astart=AST, aend=AEND):
    g = np.concatenate([
        np.logspace(np.log10(astart), 0, 20001),
        np.linspace(astart, aend, 20001),
        aend - (aend - astart) * np.logspace(-10, 0, 4001),
    ])
    return np.unique(g[(g >= astart) & (g <= aend)])


def reasons(r: dict) -> dict[str, np.ndarray]:
    """Reconstruct the physical conditions from the cached coefficients (same decision rules as the code)."""
    a = r["a"]
    out = {}
    om1 = 1 + r["EFTOmegaV"]
    beyond = (r["EFTGamma6V"] != 0) | ((r["EFTGamma3V"] + r["EFTGamma4V"]) != 0)
    if r["EFT_ghost_stability"] or r["EFT_gradient_stability"]:
        out["beyond_GLPV"] = beyond
        out["G_neg"] = om1 <= 0
    if r["EFT_ghost_stability"]:
        out["ghost_scalar"] = r["EFT_kinetic"] < 0
        out["ghost_tensor"] = r["EFTAT"] < 0
    if r["EFT_gradient_stability"]:
        out["grad_scalar"] = r["EFT_gradient"] < 0
        out["grad_tensor"] = r["EFTDT"] < 0
    if r["EFT_mass_stability"]:
        thr = -(r["EFT_mass_stability_rate"] * r["adotoa"] / a) ** 2
        out["mass"] = ((r["EFT_mu1"] < thr) | (r["EFT_mu2"] < thr)) & (a > r["pert_turn_on"])
    if r["EFT_positivity_bounds"]:
        out["posbounds"] = (r["posbound1"] < 0) | (r["posbound2"] < 0)
    if r["EFT_minkowski_limit"]:
        out["minkowski"] = (-r["ghostpos"] <= 0) | (-r["tachpos"] < 0)
    nanf = np.zeros(len(a), dtype=bool)
    for f in t6lib.DEFAULT_FIELDS:
        nanf |= ~np.isfinite(r[f])
    out["nan_mapped_fields"] = nanf
    return out


def run_one(tree, ini: Path, outdir: Path, overrides: dict) -> dict:
    ov = {"EFTCAMB_skip_stability": "T", "feedback_level": "0", **overrides}
    r0 = t6lib.probe(tree, ini, overrides=ov, a_grid=[1.0], feedback=0, name=ini.stem + "_pre")
    ast = float(r0.get("stability_time", AST))
    lin, log1, log2 = sampler_grid(ast, AEND)
    fine = fine_grid(ast, AEND)
    grid = np.unique(np.concatenate([lin, log1, fine]))
    r = t6lib.probe(tree, ini, overrides=ov, a_grid=grid, check=True, feedback=0, name=ini.stem,
                    cleanup_each=False, timeout=1800)
    summ: dict = {"model": ini.stem, "status": r["status"], "bg_error": r.get("bg_error"),
                  "full_check": r.get("full_check"), "dir": r["dir"]}
    if r["status"] != "ok" or "verdict" not in r:
        summ["stderr"] = r.get("stderr", [])[-8:]
        return summ
    a, v = r["a"], r["verdict"]
    np.savez_compressed(outdir / f"{ini.stem}.npz", **{k: r[k] for k in ["a", "verdict"] + t6lib.DEFAULT_FIELDS})
    insamp = np.isin(a, np.concatenate([lin, log1]))
    summ["pert_turn_on"] = r["pert_turn_on"]
    summ["stability_time"] = ast
    summ["flags"] = {k: r[k] for k in ("EFT_ghost_math_stability", "EFT_mass_math_stability", "EFT_ghost_stability",
                                       "EFT_gradient_stability", "EFT_mass_stability", "EFT_additional_priors",
                                       "EFT_positivity_bounds", "EFT_minkowski_limit", "EFT_mass_stability_rate")}
    summ["n_points"] = int(len(a))
    summ["n_unstable_fine"] = int((~v).sum())
    summ["n_unstable_sampler"] = int((~v[insamp]).sum())
    if (~v).any():
        bad = a[~v]
        summ["unstable_range"] = [float(bad.min()), float(bad.max())]
    rs = reasons(r)
    summ["reasons"] = {k: int(m.sum()) for k, m in rs.items() if m.any()}
    py_bad = np.zeros(len(a), dtype=bool)
    for k, m in rs.items():
        if k != "nan_mapped_fields":
            py_bad |= m
    summ["n_python_unstable"] = int(py_bad.sum())
    summ["n_fortran_only_unstable"] = int((~v & ~py_bad).sum())  # additional priors / NaN / math conditions
    summ["n_python_only_unstable"] = int((v & py_bad).sum())
    if (~v & ~py_bad).any():
        summ["fortran_only_range"] = [float(a[~v & ~py_bad].min()), float(a[~v & ~py_bad].max())]
    if (v & py_bad).any():
        summ["python_only_range"] = [float(a[v & py_bad].min()), float(a[v & py_bad].max())]
    summ["missed_by_sampler"] = bool(summ["full_check"] and (~v).any())
    return summ


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tree", default="A0")
    p.add_argument("--pattern", default="*.ini")
    p.add_argument("--inis", nargs="*", default=None)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--tag", default="scan_A0")
    p.add_argument("--set", nargs="*", default=[])
    a = p.parse_args()
    inis = [Path(x) for x in a.inis] if a.inis else harness.model_inis(a.pattern)
    outdir = t6lib.WORK / a.tag
    outdir.mkdir(parents=True, exist_ok=True)
    ov = dict(s.split("=", 1) for s in a.set)
    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(lambda ini: run_one(a.tree, ini, outdir, ov), inis))
    for s in res:
        (outdir / f"{s['model']}.json").write_text(json.dumps(s, indent=1))
        print(json.dumps({k: s.get(k) for k in ("model", "status", "full_check", "n_unstable_fine",
                                                "n_unstable_sampler", "unstable_range", "reasons",
                                                "n_fortran_only_unstable", "n_python_only_unstable",
                                                "missed_by_sampler", "bg_error")}))


if __name__ == "__main__":
    main()

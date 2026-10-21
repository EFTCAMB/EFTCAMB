"""Repro (t2d-static): full-mapping models whose background EFT functions read adotoa/Hdot/Hdotdot from the timestep
cache get the values of the PREVIOUS derivs call of that k mode (Beyond Horndeski, Extended Galileon).

Code: fortran/equations.f90:2813-2817 (derivs, full models):
        call CP%EFTCAMB%model%compute_background_EFT_functions(a, ...)   ! reads eft_cache%adotoa, Hdot, Hdotdot
        call CP%EFTCAMB%model%compute_adotoa(a, ...)                     ! only now adotoa is set for this tau
      ... compute_H_derivs (equations.f90:2861, 2893) sets Hdot, Hdotdot; compute_background_EFT_functions is called
      again only for designer models (equations.f90:2897).  EV%eft_cache is not re-initialised between derivs calls in
      release builds (only when DebugEFTCAMB, equations.f90:1955), so
  * 008p5_Beyond_Horndeski.f90:1009/1014/1051 (EFTc, EFTLambda, EFTcdot, EFTLambdadot built from adotoa, Hdot, Hdotdot)
  * 008p7_Extended_Galileon.f90:337-366 (the 'adotoa == 0' guard recomputes the background only on a fresh cache)
use the background of the previous call: the EFT functions entering the perturbation equations at time tau are not
the functions of the background at tau, and depend on the integrator's step history.
Same call order in EFTCAMB_GetEFTFunctions (11_EFTCAMB_python.f90:421-423) where the cache is zero (python-allocated).

Method: the audit tree work/t2d-static/LAG is FIX3 (= A0 + the Extended Galileon neutrino fix, which does not touch
Beyond Horndeski) plus instrumentation after the second compute_H_derivs call in derivs (full models only): it
re-evaluates compute_background_EFT_functions with the current background, prints to stderr the largest relative
difference between the value derivs used and the fresh value (a >= 1e-3), and restores the used values unless
EFT_LAG_FIX=1 (then the fresh values are kept: candidate fix).  Model: audit ini 008p5-beyond-horndeski.ini
(BH2 of arXiv:1905.11364) with EFTCAMB_skip_stability = T, 1 OpenMP thread.
Runs: A0 (unmodified), LAG (instrumented, used values kept), LAG with EFT_LAG_FIX=1.
Expected: used == fresh (relative difference 0), and the candidate fix does not change the spectra.

Usage: /Users/marco/Software/miniforge3/bin/python fullmap-bg-eft-stale-background.py
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
A0 = AUDIT / "trees" / "A0"
LAG = AUDIT / "work" / "t2d-static" / "LAG"
INI = Path(__file__).resolve().parents[1] / "inis" / "008p5-beyond-horndeski.ini"
FIELDS = ["EFTc", "EFTLambda", "EFTcdot", "EFTLambdadot"]

CHILD = r"""
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
p = camb.read_ini(sys.argv[2], no_validate=True)
r = camb.get_results(p)
cl = r.get_cmb_power_spectra(p, CMB_unit="muK", raw_cl=True)
k, z, pk = r.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
np.savez(sys.argv[3], tt=cl["unlensed_scalar"][2:2501, 0], ee=cl["unlensed_scalar"][2:2501, 1],
         pp=cl["lens_potential"][2:2501, 0], k=k, pk=pk[0])
"""


def run(tree: Path, tmp: Path, name: str, env_extra: dict[str, str]) -> tuple[dict, str]:
    w = tmp / f"{name}.ini"
    w.write_text(f"DEFAULT({INI})\nEFTCAMB_skip_stability = T\nEFTCAMB_write_background = F\n"
                 f"feedback_level = 0\nEFTCAMB_feedback_level = 0\noutput_root = {tmp / name}\n")
    out = tmp / f"{name}.npz"
    env = dict(os.environ, OMP_NUM_THREADS="1", **env_extra)
    p = subprocess.run([sys.executable, "-c", CHILD, str(tree), str(w), str(out)], capture_output=True, text=True,
                       env=env, timeout=900, cwd=tmp)
    if p.returncode:
        raise RuntimeError(f"{name}: rc={p.returncode}\n{p.stderr[-2000:]}")
    import numpy as np

    return dict(np.load(out)), p.stderr


def lagmax(stderr: str) -> dict[str, tuple[float, float, float, float, float]]:
    best: dict[str, tuple] = {}
    for line in stderr.splitlines():
        if line.startswith("LAGMAX"):
            v = line.split()
            i, (a, k, tau, used, fresh, rel) = int(v[1]), map(float, v[2:8])
            f = FIELDS[i - 1]
            if f not in best or rel > best[f][-1]:
                best[f] = (a, k, used, fresh, rel)
    return best


def maxrel(x, y) -> float:
    import numpy as np

    return float(np.max(np.abs(x / y - 1)))


def main() -> None:
    try:  # share the audit's machine-wide CPU slots when the harness is available
        sys.path.insert(1, str(Path(__file__).resolve().parents[1] / "checks"))
        from harness import slot
    except Exception:  # noqa: BLE001
        import contextlib

        def slot(n):  # noqa: ARG001
            return contextlib.nullcontext()

    from concurrent.futures import ThreadPoolExecutor

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)

        def job(args):
            tree, name, env = args
            with slot(1):
                return run(tree, tmp, name, env)

        with ThreadPoolExecutor(3) as ex:
            (r_a0, _), (r_lag, err_lag), (r_fix, err_fix) = ex.map(
                job, [(A0, "a0", {}), (LAG, "lag", {}), (LAG, "fix", {"EFT_LAG_FIX": "1"})])
    print(f"model: {INI.name} (Beyond Horndeski, EFTflag=4, FullMappingEFTmodel=5), EFTCAMB_skip_stability = T")
    print("control, A0 vs LAG with the used values kept (only recompilation/rounding differs):")
    print(f"  unlensed TT max rel diff = {maxrel(r_a0['tt'], r_lag['tt']):.3e}   EE = {maxrel(r_a0['ee'], r_lag['ee']):.3e}"
          f"   C_l^phiphi = {maxrel(r_a0['pp'], r_lag['pp']):.3e}   P(k, z=0) = {maxrel(r_a0['pk'], r_lag['pk']):.3e}")
    print("EFT functions used by derivs (background of the previous call) vs re-evaluated with the current background:")
    for f, (a, k, used, fresh, rel) in lagmax(err_lag).items():
        print(f"  {f:13s} max rel diff = {rel:.3e}   at a = {a:.4f}, k = {k:.3e}/Mpc: used {used:+.5e}, fresh {fresh:+.5e}")
    print("candidate fix (fresh values kept, EFT_LAG_FIX=1) vs code under audit:")
    print(f"  unlensed TT max rel diff = {maxrel(r_fix['tt'], r_lag['tt']):.3e}   EE = {maxrel(r_fix['ee'], r_lag['ee']):.3e}"
          f"   C_l^phiphi = {maxrel(r_fix['pp'], r_lag['pp']):.3e}   P(k, z=0) = {maxrel(r_fix['pk'], r_lag['pk']):.3e}")
    print("expected: max rel diff used/fresh = 0 for every field and no change of the spectra with the candidate fix")


if __name__ == "__main__":
    main()

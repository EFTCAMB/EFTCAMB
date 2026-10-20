"""Repro: reusing the parameters of an EFTCAMB result (results.Params) aborts the process.

fortran/results.f90:580 (CAMBdata_SetParams, EFTCAMB block):
    call this%CP%eft_par_cache%initialize()          ! does not deallocate grhormass / nu_masses
    ...
    allocate(this%CP%eft_par_cache%grhormass(max_nu), this%CP%eft_par_cache%nu_masses(max_nu))
After a calculation with EFTflag /= 0, State%CP%eft_par_cache%grhormass is allocated.  Any later SetParams whose
input CAMBparams is a copy of State%CP (this%CP = P copies the allocated arrays) hits the allocate statement with
allocated arrays and gfortran stops the whole Python process ("Attempting to allocate already allocated variable").

Reachable through documented CAMB v2 APIs, e.g. results.get_BAO(z) (params defaults to self.Params.copy()),
camb.get_background(results.Params), results.calc_background(results.Params.copy()).
GR (EFTflag = 0) is unaffected.

Expected: the second call returns (same H(z=1)); observed: Fortran runtime error, exit code 2.
Run: /Users/marco/Software/miniforge3/bin/python eft-par-cache-realloc-abort.py [tree]
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)


def child(tree: str, eft: bool, mode: str) -> None:
    sys.path.insert(0, tree)
    import camb

    kw = dict(H0=70, ombh2=0.0226, omch2=0.112, mnu=0.06)
    if eft:
        kw.update(EFTflag=1, PureEFTmodel=1, PureEFTmodelOmega=1, EFTwDE=0, EFTOmega0=0.1, feedback_level=0)
    r = camb.get_background(camb.set_params(**kw))
    print(f"    first calculation: H(z=1) = {r.hubble_parameter(1.0):.12f}", flush=True)
    if mode == "get_BAO":
        bao = r.get_BAO([0.5, 1.0])
        print(f"    results.get_BAO([0.5, 1]): H = {bao[:, 1]}", flush=True)
    else:
        r2 = camb.get_background(r.Params)
        print(f"    camb.get_background(results.Params): H(z=1) = {r2.hubble_parameter(1.0):.12f}", flush=True)


def main() -> None:
    arg = sys.argv[1] if len(sys.argv) > 1 else "A0"
    tree = arg if os.path.isabs(arg) else str(ROOT / "trees" / arg)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    print(f"tree: {tree}")
    failures = 0
    for eft in (False, True):
        for mode in ("get_BAO", "reuse_Params"):
            p = subprocess.run([sys.executable, __file__, "_child", tree, "1" if eft else "0", mode], env=env,
                               capture_output=True, text=True)
            print(f"  EFTflag={'1' if eft else '0'} {mode}: exit code {p.returncode}")
            print(p.stdout.rstrip())
            err = [ln for ln in p.stderr.splitlines() if "runtime error" in ln or ln.startswith("At line")]
            for ln in err:
                print("    " + ln)
            failures += p.returncode != 0
    print("expected: all exit code 0; observed:", "OK" if failures == 0 else f"{failures} aborted")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_child":
        child(sys.argv[2], sys.argv[3] == "1", sys.argv[4])
    else:
        main()

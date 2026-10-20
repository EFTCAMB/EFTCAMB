"""Repro: EV%EFTCAMB_printing (and EV%EFTCAMB_QS) are garbage in the OpenMP source loop, so a random subset of k modes
runs every derivs call in "printing" mode.

Code: equations.f90:195-196 give the EvolutionVars components default initialisers (EFTCAMB_printing = .false.,
EFTCAMB_QS = .false.), but cmbmain.f90:249 evolves the k modes with !$OMP PARALLEL DO ... PRIVATE(EV): gfortran does
NOT default-initialise PRIVATE copies of derived types, and nothing in DoSourcek/initial()/GetNumEqns resets these two
flags (upstream CAMB nullifies its own pointer components in initial(), equations.f90:2140-2142, for this reason).
With EFTCAMB_printing = .true. every derivs call also runs compute_tensor_factors, compute_stability_factors,
compute_mass_factors (equations.f90:3628-3634; mass factors read never-computed EFTcdotdot/EFTGamma*PP), the whole
output block incl. the QS mu/Sigma formulas (equations.f90:3638-3768) and the cache copies.

Method:
 (1) a 20-line Fortran program, compiled with the same gfortran/-fopenmp, shows that a PRIVATE copy of a derived type
     with default initialisation starts with garbage (and keeps it between iterations);
 (2) tree work/t1-crash/relP = A0 (release) + one print in DoSourcek of the incoming EV%EFTCAMB_printing:
     count of k modes that start with EFTCAMB_printing = .true.;
 (3) tree work/t1-crash/relQ = A0 + "EV%EFTCAMB_printing = .false.; EV%EFTCAMB_QS = .false." at the top of DoSourcek:
     CPU time and spectra compared with A0 (same model, same threads).
Expected: (1) F, (2) 0 modes, (3) identical cost.

Run: /Users/marco/Software/miniforge3/bin/python omp-private-ev-printing-uninit.py
"""

import json
import os
import resource
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

try:  # hold a machine-wide CPU slot when the audit harness is available (shared machine); optional otherwise
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
    import harness

    def _slot(n):
        return harness.slot(n)
except Exception:  # noqa: BLE001
    from contextlib import nullcontext

    def _slot(n):
        return nullcontext()

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
             "scratchpad/audit")
W = AUDIT / "work" / "t1-crash"
PAR = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")
GFORTRAN = "/Users/marco/Software/miniforge3/bin/gfortran"

FSRC = """
module m
  implicit none
  type ev_t
    real(8) :: q
    logical :: printing = .false.
  end type
end module
program p
  use m
  implicit none
  type(ev_t) :: ev
  integer :: i
  !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(DYNAMIC), PRIVATE(ev, i)
  do i = 4, 1, -1
    print '(a,i2,a,l2)', ' iteration', i, ': private ev%printing =', ev%printing
  end do
  !$OMP END PARALLEL DO
end program
"""

CHILD = r"""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_results(cp)
cl = res.get_cmb_power_spectra(cp, CMB_unit='muK', raw_cl=True)['total']
np.save(sys.argv[3], cl)
"""


def fortran_demo() -> None:
    with tempfile.TemporaryDirectory() as d:
        src = Path(d) / "t.f90"
        src.write_text(FSRC)
        r = subprocess.run([GFORTRAN, "-fopenmp", "-O0", str(src), "-o", str(Path(d) / "t")], capture_output=True,
                           text=True)
        if r.returncode:
            print("gfortran failed:", r.stderr[-300:])
            return
        out = subprocess.run([str(Path(d) / "t")], env=dict(os.environ, OMP_NUM_THREADS="1"), capture_output=True,
                             text=True).stdout
        print("(1) gfortran PRIVATE copy of a default-initialised component (expected F every iteration):")
        print("\n".join("   " + ln.strip() for ln in out.strip().splitlines()))


def run_child(tree: Path, model: str, out: Path, threads: int) -> tuple[float, str, np.ndarray]:
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    env.pop("PYTHONPATH", None)
    with _slot(threads):
        before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
        r = subprocess.run([sys.executable, "-c", CHILD, str(tree), str(PAR / f"{model}.ini"), str(out)], env=env,
                           capture_output=True, text=True, cwd=out.parent, timeout=3600, preexec_fn=lambda: os.nice(10))
        cpu = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime - before
    if r.returncode:
        raise SystemExit(f"run failed: {r.stderr[-500:]}")
    return cpu, r.stdout, np.load(out)


def main() -> None:
    fortran_demo()
    models = sys.argv[1:] or ["2_PEFT_Omega_const_1"]
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        print("(2) k modes entering DoSourcek with EV%EFTCAMB_printing = .true. (release build + probe print):")
        for m in models:
            _, out, _ = run_child(W / "relP" / "A0", m, d / f"p_{m}.npy", 2)
            print(f"   {m:24s}: {out.count('AUDIT_PRINTING_TRUE'):4d} true, {out.count('AUDIT_PRINTING_FALSE'):4d} false"
                  "  (expected 0 true)")
        print("(3) A0 (code under audit) vs relQ (flags reset in DoSourcek): CPU seconds, 1 thread, total Cls")
        for m in models:
            ta, _, ca = run_child(AUDIT / "trees" / "A0", m, d / f"a_{m}.npy", 1)
            tq, _, cq = run_child(W / "relQ" / "A0", m, d / f"q_{m}.npy", 1)
            n = min(len(ca), len(cq), 2501)
            rel = np.max(np.abs(ca[2:n, :3] - cq[2:n, :3]) / np.maximum(np.abs(cq[2:n, :3]), 1e-300))
            print(f"   {m:24s}: CPU A0 = {ta:7.1f} s, fixed = {tq:7.1f} s (A0/fixed = {ta / tq:5.2f}); "
                  f"max |dCl/Cl| TT/EE/BB = {rel:.2e}, bit-identical = {np.array_equal(ca, cq)}")


if __name__ == "__main__":
    main()

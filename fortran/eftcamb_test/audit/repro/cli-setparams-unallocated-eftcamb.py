"""Repro: the command-line camb executable segfaults in CAMBdata_SetParams (results.f90:365) for every model when built
without short-circuit evaluation (gfortran -O0 / -g debug builds).

Code (EFTCAMB MOD block, results.f90:364-368):
    if (.not. allocated(P%DarkEnergy) .and. this%CP%EFTCAMB%EFTflag == 0) then
this%CP is the State's *previous* parameter set; it is assigned from P only ~20 lines later (this%CP = P), so on the
first call from the command line this%CP%EFTCAMB is an unallocated allocatable.  Fortran .and. does not
short-circuit: gfortran -O0 evaluates both operands and dereferences address 0 (EXC_BAD_ACCESS at results.f90:365).
Release builds survive only because the optimiser happens to skip the second operand (P%DarkEnergy is allocated).
The test is also logically wrong: it should look at the new parameters P, not this%CP.

Method: run 1_EFT_GR.ini with
  (a) the command-line exe of the code under audit built with the AD flags (-O0 -fcheck=all; plus only the F-002
      upstream patch), tree work/t1-crash/chk/AD;
  (b) the same flags plus the proposed fix (test P%EFTCAMB, guarded by allocated()), tree work/t1-crash/chk3x/AD;
  (c) the release exe of the code under audit (tree work/t1-crash/rexe/A0).
Expected: all three start the calculation.  Observed: (a) SIGSEGV within a second.

Run: /Users/marco/Software/miniforge3/bin/python cli-setparams-unallocated-eftcamb.py
"""

import os
import shutil
import signal
import subprocess
import tempfile
import time

W = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/t1-crash"
INI = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/1_EFT_GR.ini"
EXES = {
    "(a) -O0 checked exe, code under audit": f"{W}/chk/AD/fortran/camb",
    "(b) -O0 checked exe + fix": f"{W}/chk3x/AD/fortran/camb",
    "(c) release exe, code under audit": f"{W}/rexe/A0/fortran/camb",
}
RUN_SECONDS = 60  # long enough to get past SetParams; the run is then stopped


def main() -> None:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    for label, exe in EXES.items():
        if not os.path.exists(exe):
            print(f"{label}: {exe} not built, skipped")
            continue
        with tempfile.TemporaryDirectory() as d:
            shutil.copy(os.path.join(os.path.dirname(exe), "HighLExtrapTemplate_lenspotentialCls.dat"), d)
            with open(os.path.join(d, "run.ini"), "w") as f:
                f.write(f"DEFAULT({INI})\noutput_root = {d}/gr\nnumber_of_threads = 1\nget_tensor_cls = F\n")
            t0 = time.time()
            p = subprocess.Popen([exe, "run.ini"], cwd=d, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                 text=True)
            try:
                out, _ = p.communicate(timeout=RUN_SECONDS)
                rc = p.returncode
            except subprocess.TimeoutExpired:
                p.kill()
                out, _ = p.communicate()
                rc = None
            dt = time.time() - t0
        if rc is None:
            obs = f"still running after {RUN_SECONDS}s (past SetParams) -> ok"
        elif rc < 0:
            obs = f"killed by {signal.Signals(-rc).name} after {dt:.1f}s"
        else:
            obs = f"exit code {rc} after {dt:.1f}s"
        sig = [ln for ln in out.splitlines() if "Program received signal" in ln]
        print(f"{label:42s}: expected ok, observed {obs} {sig[0] if sig else ''}")


if __name__ == "__main__":
    main()

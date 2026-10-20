"""Build the command-line camb with the *shipped* gfortran Debug target (fortran/Makefile DEBUGFLAGS:
-g -fbacktrace -ffpe-trap=invalid,overflow,zero -fbounds-check -DDEBUG, default -O0) from a copy of the frozen A0 tree.

Variants:
  dbg     : unmodified A0 sources (no F-002 patch: the SetParams crash, if any, happens before thermo init)
  dbgfix  : A0 + only the proposed fix of results.f90:365 (nothing else changed)
  dbgprobe: A0 + prints of the allocation status just before / after the check (check unchanged)

Usage: python cli_setparams_build_debug.py dbg|dbgfix|dbgprobe
Holds the machine-wide build lock and 4 CPU slots, nice 10, make -j4.
"""

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

S = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45"
    "/scratchpad/audit"
)
OUT = S / "work" / "verify" / "cli_setparams"
DEBUGFLAGS = ("-g -fbacktrace -ffpe-trap=invalid,overflow,zero -fbounds-check -MMD -cpp -ffree-line-length-none "
              "-fallow-argument-mismatch -fmax-errors=4 -fopenmp -lblas -llapack -DDEBUG")

OLD = """    if (.not. allocated(P%DarkEnergy) .and. this%CP%EFTCAMB%EFTflag == 0) then
        call GlobalError('DarkEnergy not set', error_darkenergy)
    end if
"""
NEW = """    if (.not. allocated(P%DarkEnergy)) then
        if (.not. allocated(P%EFTCAMB)) then
            call GlobalError('DarkEnergy not set', error_darkenergy)
        else if (P%EFTCAMB%EFTflag == 0) then
            call GlobalError('DarkEnergy not set', error_darkenergy)
        end if
    end if
"""

PROBE_PRE = """    write(*, '(a,3l2)') 'AUDITPROBE allocated(P%DarkEnergy), allocated(this%CP%EFTCAMB), allocated(P%EFTCAMB):', &
        allocated(P%DarkEnergy), allocated(this%CP%EFTCAMB), allocated(P%EFTCAMB)
    flush(6)
"""
PROBE_POST = """    write(*, '(a)') 'AUDITPROBE passed the DarkEnergy check'
    flush(6)
"""


def main(variant: str) -> None:
    dest = OUT / variant
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(S / "trees" / "A0", dest, symlinks=True)
    exe = dest / "fortran" / "camb"
    if exe.exists():
        exe.unlink()
    res = dest / "fortran" / "results.f90"
    txt = res.read_text()
    assert txt.count(OLD) == 1, "pattern not found exactly once"
    if variant == "dbgfix":
        res.write_text(txt.replace(OLD, NEW))
    elif variant == "dbgprobe":
        res.write_text(txt.replace(OLD, PROBE_PRE + OLD + PROBE_POST))
    elif variant != "dbg":
        raise SystemExit("variant must be dbg, dbgfix or dbgprobe")
    build = Path(f"{dest}_build")
    if build.exists():
        shutil.rmtree(build)
    env = dict(os.environ, CAMB_BUILD_DIR=str(build), MAKEFLAGS="-j4")
    log = Path(f"{dest}_build.log")
    t0 = time.time()
    # forutils' Makefile has no module dependencies, so its parallel build can race: build forutils serially
    # first (-j1, attempt 1), then the rest with -j4 (attempts 2-3 resume the same build directory).
    with harness._build_lock(), harness.slot(4), open(log, "w") as f:
        for attempt, jobs in enumerate(["-j1", "-j4", "-j4"], 1):
            env["MAKEFLAGS"] = jobs
            targets = ["Debug"] if attempt > 1 else ["Debug", "libforutils", "libforutils_so"]
            if attempt == 1:
                targets = ["libforutils", "libforutils_so", "OUTPUT_DIR=" + str(build / "camb" / "Debug"),
                           "DLL_DIR=" + str(build / "camb" / "Debuglib"), "FFLAGS=" + DEBUGFLAGS]
            f.write(f"\n===== attempt {attempt}: make {' '.join(targets)} ({jobs})\n")
            f.flush()
            rc = subprocess.run(["make", *targets], cwd=dest / "fortran", env=env, stdout=f,
                                stderr=subprocess.STDOUT, preexec_fn=lambda: os.nice(10)).returncode
            if attempt > 1 and rc == 0:
                break
    print(f"variant={variant} rc={rc} exe_exists={exe.exists()} elapsed={time.time() - t0:.1f}s log={log}")
    # record the flags actually used for results.f90
    for ln in log.read_text().splitlines():
        if "results.f90" in ln and "-c" in ln:
            print("compile line:", ln.strip()[:400])
            break


if __name__ == "__main__":
    main(sys.argv[1])

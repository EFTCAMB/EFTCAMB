"""Independent check of finding cli-setparams-unallocated-eftcamb.

Uses the *shipped* gfortran `make Debug` configuration (fortran/Makefile DEBUGFLAGS, -DDEBUG, default -O0, no
-fcheck=all, no F-002 patch) built by cli_setparams_build_debug.py from a copy of the frozen A0 tree:
  dbg    : unmodified code under audit
  dbgfix : only the proposed results.f90:365 fix applied
For several models (EFTflag 0, pure EFT, designer, full mapping, and the upstream CAMB inifiles/params.ini without
any EFTCAMB key) it runs the command-line exe and records how it ends, and runs lldb with a breakpoint on
CheckSupportedParams (the statement right after the check in CAMBdata_SetParams): reaching it means SetParams
passed the check.

Usage: python cli_setparams_independent.py [runs|probe|lldb|all]
"""

import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

REPO = Path("/Users/marco/IW/eftcamb_update/EFTCAMB")
V = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45"
    "/scratchpad/audit/work/verify/cli_setparams"
)
PAR = REPO / "fortran" / "eftcamb_test" / "parameters"
MODELS = {
    "1_EFT_GR(EFTflag=0)": PAR / "1_EFT_GR.ini",
    "2_PEFT_Omega_const_1": PAR / "2_PEFT_Omega_const_1.ini",
    "4_Des_F_R_LCDM_1": PAR / "4_Des_F_R_LCDM_1.ini",
    "5_quint_exponential": PAR / "5_quint_exponential.ini",
    "upstream inifiles/params.ini": REPO / "inifiles" / "params.ini",
}
RUN_SECONDS = 45


def prep(dirname: Path, ini: Path, exe: Path) -> Path:
    if dirname.exists():
        shutil.rmtree(dirname)
    dirname.mkdir(parents=True)
    shutil.copy(exe.parent / "HighLExtrapTemplate_lenspotentialCls.dat", dirname)
    (dirname / "run.ini").write_text(
        f"DEFAULT({ini})\noutput_root = {dirname}/out\nnumber_of_threads = 1\nget_tensor_cls = F\n"
        "EFTCAMB_write_background = F\n"
    )
    return dirname / "run.ini"


def run_group(cmd: list[str], cwd: Path, timeout: float) -> tuple[int | None, str, float]:
    t0 = time.time()
    p = subprocess.Popen(cmd, cwd=cwd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, start_new_session=True, env=dict(os.environ, OMP_NUM_THREADS="1"),
                         preexec_fn=lambda: os.nice(10))
    try:
        out, _ = p.communicate(timeout=timeout)
        rc = p.returncode
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        out, _ = p.communicate()
        rc = None
    return rc, out, time.time() - t0


def describe(rc: int | None, dt: float) -> str:
    if rc is None:
        return f"still running after {dt:.0f}s (killed)"
    if rc < 0:
        return f"killed by {signal.Signals(-rc).name} after {dt:.1f}s"
    return f"exit code {rc} after {dt:.1f}s"


def runs() -> None:
    print("=== plain command-line runs (shipped make Debug flags) ===")
    for variant in ("dbg", "dbgfix"):
        exe = V / variant / "fortran" / "camb"
        for label, ini in MODELS.items():
            d = V / "runs" / variant / re.sub(r"[^A-Za-z0-9_]", "_", label)
            prep(d, ini, exe)
            with harness.slot(1):
                rc, out, dt = run_group([str(exe), "run.ini"], d, RUN_SECONDS)
            (d / "stdout.txt").write_text(out)
            sig = [ln.strip() for ln in out.splitlines() if "Program received signal" in ln]
            frame = [ln.strip() for ln in out.splitlines() if re.search(r"results\.f90:\d+", ln)]
            print(f"{variant:7s} {label:30s}: {describe(rc, dt):34s} {sig[0] if sig else ''} "
                  f"{frame[0] if frame else ''}")


LLDB_CMDS = [
    "br set -n __results_MOD_checksupportedparams",
    "run",
    "bt 3",
    "frame select 0",
    "p this->cp.eftcamb._data",
    "p p->darkenergy._data",
    "p p->eftcamb._data",
    "kill",
]


def lldb() -> None:
    print("=== lldb: breakpoint on CheckSupportedParams (next statement after the check) ===")
    for variant in ("dbg", "dbgfix"):
        exe = V / variant / "fortran" / "camb"
        for label in ("1_EFT_GR(EFTflag=0)", "4_Des_F_R_LCDM_1"):
            d = V / "lldb" / variant / re.sub(r"[^A-Za-z0-9_]", "_", label)
            prep(d, MODELS[label], exe)
            cmd = ["lldb", "--batch"]
            for c in LLDB_CMDS:
                cmd += ["-o", c]
            cmd += ["--", str(exe), "run.ini"]
            with harness.slot(1):
                rc, out, dt = run_group(cmd, d, 240)
            (d / "lldb.txt").write_text(out)
            keep = [ln for ln in out.splitlines() if re.search(
                r"stop reason|frame #|EXC_BAD|results\.f90|\(.*\*\) \$|= 0x|error:", ln)]
            print(f"--- {variant} {label}: lldb {describe(rc, dt)}")
            for ln in keep[:14]:
                print("   ", ln.strip()[:220])


def probe() -> None:
    print("=== dbgprobe: allocation status printed just before the unchanged check, marker printed after it ===")
    exe = V / "dbgprobe" / "fortran" / "camb"
    for label in ("1_EFT_GR(EFTflag=0)", "4_Des_F_R_LCDM_1", "upstream inifiles/params.ini"):
        d = V / "runs" / "dbgprobe" / re.sub(r"[^A-Za-z0-9_]", "_", label)
        prep(d, MODELS[label], exe)
        with harness.slot(1):
            rc, out, dt = run_group([str(exe), "run.ini"], d, RUN_SECONDS)
        (d / "stdout.txt").write_text(out)
        keep = [ln.strip() for ln in out.splitlines() if "AUDITPROBE" in ln or "Program received" in ln]
        print(f"{label:30s}: {describe(rc, dt)}")
        for ln in keep:
            print("    ", ln)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("runs", "all"):
        runs()
    if what == "lldb":  # needs debugger attach rights (macOS Developer mode); not part of "all"
        lldb()
    if what in ("probe", "all"):
        probe()

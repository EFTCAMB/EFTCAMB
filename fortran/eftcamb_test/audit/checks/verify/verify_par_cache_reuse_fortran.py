"""Independent pure-Fortran check of finding eft-par-cache-realloc-abort.

Compiles par_cache_reuse_ftest.f90 against the frozen build objects of a tree (read-only use of the .mod/.o files)
and calls CAMBdata%SetParams twice: first with user parameters read by CAMB_ReadParamFile, then with
P2 = State%CP (mode reuse), with P2 after deallocating grhormass/nu_masses (reuse_dealloc), or with the original P
(user). No Python wrapper is involved.

Usage: python verify_par_cache_reuse_fortran.py [build_label ...]   (labels: A0, fix)
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
HERE = Path(__file__).resolve().parent
PARAMS = HERE.parents[2] / "parameters"
OUT = ROOT / "work" / "verify" / "par_cache_realloc" / "ftest"
BUILDS = {
    "A0": ROOT / "build" / "A0",
    "fix": ROOT / "work" / "verify" / "fixmin_realloc_build",
}
MODELS = {
    "fR": ("4_Des_F_R_LCDM_1", ""),
    "PEFT": ("2_PEFT_Omega_const_1", ""),
    "GR": ("4_Des_F_R_LCDM_1", "EFTflag = 0\n"),
}


def compile_test(label: str) -> Path:
    b = BUILDS[label]
    lib, fu = b / "camb" / "Releaselib", b / "forutils" / "Releaselib"
    exe = OUT / f"ftest_{label}"
    objs = sorted(str(p) for p in lib.glob("*.o"))
    cmd = ["gfortran", "-O1", "-cpp", "-ffree-line-length-none", "-fopenmp", "-fbacktrace", f"-I{lib}", f"-I{fu}",
           "-J", str(OUT / f"mod_{label}"), str(HERE / "par_cache_reuse_ftest.f90"), *objs,
           str(fu / "libforutils.a"), "-lblas", "-llapack", "-o", str(exe)]
    (OUT / f"mod_{label}").mkdir(parents=True, exist_ok=True)
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode:
        print(p.stderr[-3000:])
        raise SystemExit(f"compile failed for {label}")
    return exe


def make_ini(name: str) -> Path:
    model, extra = MODELS[name]
    ini = OUT / f"{name}.ini"
    ini.write_text(f"DEFAULT({PARAMS / (model + '.ini')})\noutput_root = {OUT}/{name}\n"
                   f"EFTCAMB_write_background = F\nfeedback_level = 0\nEFTCAMB_feedback_level = 0\n{extra}")
    return ini


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    labels = sys.argv[1:] or ["A0", "fix"]
    env = dict(os.environ, OMP_NUM_THREADS="1")
    for label in labels:
        exe = compile_test(label)
        for name in MODELS:
            ini = make_ini(name)
            for mode in ("reuse", "reuse_dealloc", "user"):
                with harness.slot(1):
                    p = subprocess.run(["nice", "-n", "10", str(exe), str(ini), mode], capture_output=True, text=True,
                                       env=env, cwd=str(OUT), timeout=600)
                lines = [ln for ln in (p.stdout + p.stderr).splitlines() if ln.strip()]
                keep = [ln for ln in lines if any(s in ln for s in ("EFTflag", "allocated", "H(z=1)", "error",
                                                                  "At line", "DONE", "deallocated", "original"))]
                print(f"== build={label} model={name} mode={mode} rc={p.returncode}")
                for ln in keep:
                    print("   " + ln.strip().replace(str(ROOT), "$S"))
                sys.stdout.flush()


if __name__ == "__main__":
    main()

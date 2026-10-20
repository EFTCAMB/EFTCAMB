"""Verifier for finding stability-cache-cdotdot-uninitialised (independent of the ctypes repro).

Builds cdotdot_stack.f90 against a frozen release library (modules from build/<T>/camb/Releaselib) and runs
  1. init      : which TEFTCAMB_timestep_cache fields survive %initialize() when the cache starts as NaN / 1.0
  2. direct    : EFTCAMB_Stability_Check called right after 1 MB of stack was left holding NaN vs 0 by a previous call
  3. setparams : the normal CAMBdata%SetParams path (stability forced on) after the same stack dirtying
usage: python cdotdot_indep.py [TREE=A0] [models...]
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
PAR = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")
TREE = sys.argv[1] if len(sys.argv) > 1 else "A0"
MODELS = sys.argv[2:] or ["5_quint_exponential", "4_Des_F_R_LCDM_1", "5_Kmouflage_1", "2_PEFT_gamma1_const_1",
                          "1_EFT_GR"]
W = AUDIT / "work" / "verify" / "cdotdot" / TREE
LIB = "/Users/marco/Software/miniforge3/lib"


def build() -> Path:
    W.mkdir(parents=True, exist_ok=True)
    B = AUDIT / "build" / TREE
    so = B / "camb" / "Releaselib" / "camblib.so"
    exe = W / "cdotdot_stack"
    cmd = ["gfortran", "-O0", "-g", "-cpp", "-ffree-line-length-none", "-fopenmp", "-J" + str(W), "-I" + str(W),
           f"-I{B}/camb/Releaselib", f"-I{B}/forutils/Releaselib", str(HERE / "cdotdot_stack.f90"), str(so),
           f"-Wl,-rpath,{LIB}", f"-Wl,-rpath,{so.parent}", "-o", str(exe)]
    p = subprocess.run(["nice", "-n", "10", *cmd], capture_output=True, text=True)
    if p.returncode:
        sys.exit(p.stderr[-3000:])
    return exe


def run(exe: Path, test: str, fill: str, model: str | None = None) -> tuple[list[str], list[str]]:
    with tempfile.TemporaryDirectory() as d:
        args = [str(exe), test, fill]
        if model:
            ini = Path(d) / "m.ini"
            ini.write_text(f"DEFAULT({PAR}/{model}.ini)\nfeedback_level = 3\nEFTCAMB_feedback_level = 3\n"
                           f"EFTCAMB_write_background = F\noutput_root = {d}/out\nEFTCAMB_skip_stability = F\n")
            args.append(str(ini))
        env = dict(os.environ, OMP_NUM_THREADS="1")
        with harness.slot(1):
            p = subprocess.run(["nice", "-n", "10", *args], capture_output=True, text=True, cwd=d, env=env,
                               timeout=600)
    out = p.stdout.splitlines()
    res = [l for l in out if l.startswith("RESULT")]
    nans = [l.strip() for l in out if l.strip().endswith("is NaN")]
    if not res:
        res = [f"NO RESULT rc={p.returncode} {p.stdout[-400:]} {p.stderr[-400:]}"]
    return res, nans


def main():
    exe = build()
    print("tree:", TREE, "exe:", exe)
    for fill in ("zero", "one", "nan"):
        res, nans = run(exe, "init", fill)
        print(f"[init fill={fill}]", res, "first NaN field (is_nan prints cumulatively after it):", nans[:1])
    for model in MODELS:
        for test in ("direct", "setparams"):
            row = []
            if test == "direct" and model == "1_EFT_GR":
                continue  # GR has no EFT model object; results.f90 never calls the stability check for EFTflag=0
            for fill in ("zero", "nan"):
                res, nans = run(exe, test, fill, model)
                row.append(f"{fill}: {res[0].replace('RESULT ', '')} first_nan={nans[:1]}")
            if row:
                print(f"[{test:9s}] {model:24s} " + " | ".join(row))


if __name__ == "__main__":
    main()

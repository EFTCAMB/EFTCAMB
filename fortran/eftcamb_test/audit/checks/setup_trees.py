"""Create and build frozen source trees used as audit baselines.

Each tree is a ``git archive`` export of a branch plus the forutils submodule at the commit that branch pins,
so the audit never rebuilds the main checkout.  Builds use their own CAMB_BUILD_DIR.

Trees:
    A0   eftcamb-audit start commit, release flags (the code under audit)
    AD   same commit, gfortran runtime checks (-fcheck=all, FPE traps, SNaN init), no -DDEBUG
    ADD  same commit, debug flags with -DDEBUG (EFTCAMB debug dump path)
    B0   plain CAMB v2 (branch camb)
    B1   pre-rebase EFTCAMB (branch eftcamb-pre-rebase)

Usage: python setup_trees.py ROOT [--only A0,B0] [--rebuild]
"""

from __future__ import annotations

import argparse
import io
import json
import os
import subprocess
import sys
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]

COMMON = "-MMD -cpp -ffree-line-length-none -fallow-argument-mismatch -fmax-errors=4 -fopenmp -lblas -llapack"
CHECK_FLAGS = (
    "-g -O0 -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan "
    "-finit-integer=-999999 -finit-logical=true"
)

TREES = {
    "A0": {"ref": "eftcamb-audit", "fflags": None},
    "AD": {"ref": "eftcamb-audit", "fflags": f"{CHECK_FLAGS} {COMMON}"},
    "ADD": {"ref": "eftcamb-audit", "fflags": f"{CHECK_FLAGS} {COMMON} -DDEBUG"},
    "B0": {"ref": "camb", "fflags": None},
    # the pre-rebase Makefile requests static runtimes whenever libgfortran.a exists, but conda gcc has no libgomp.a
    "B1": {"ref": "eftcamb-pre-rebase", "fflags": None, "make_args": ["SFFLAGS=-dynamiclib"]},
}


def git(*args: str, cwd: Path = REPO) -> bytes:
    return subprocess.run(["git", "-C", str(cwd), *args], check=True, capture_output=True).stdout


def export(ref: str, dest: Path) -> dict:
    commit = git("rev-parse", ref).decode().strip()
    dest.mkdir(parents=True, exist_ok=True)
    tarfile.open(fileobj=io.BytesIO(git("archive", commit)), mode="r:").extractall(dest, filter="tar")
    fu = git("ls-tree", commit, "forutils").decode().split()[2]
    tarfile.open(fileobj=io.BytesIO(git("archive", fu, cwd=REPO / "forutils")), mode="r:").extractall(
        dest / "forutils", filter="tar"
    )
    return {"ref": ref, "commit": commit, "forutils": fu}


def build(name: str, root: Path, rebuild: bool) -> tuple[str, int]:
    spec = TREES[name]
    tree = root / "trees" / name
    log = root / "logs" / f"build_{name}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    if rebuild or not tree.exists():
        info = export(spec["ref"], tree)
        info["fflags"] = spec["fflags"]
        (tree / "AUDIT_TREE.json").write_text(json.dumps(info, indent=2))
    env = dict(os.environ, CAMB_BUILD_DIR=str(root / "build" / name))
    if spec["fflags"] is None and not spec.get("make_args"):
        cmd = [sys.executable, "setup.py", "make"]
        cwd = tree
    else:
        cmd = ["make", "python", "PYCAMB_OUTPUT_DIR=../camb/", *spec.get("make_args", [])]
        if spec["fflags"] is not None:
            cmd.append(f"FFLAGS={spec['fflags']}")
        cwd = tree / "fortran"
    with open(log, "w") as f:
        rc = subprocess.run(
            cmd, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT, preexec_fn=lambda: os.nice(10)
        ).returncode
    ok = rc == 0 and (tree / "camb" / "camblib.so").exists()
    return name, 0 if ok else (rc or 1)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("root", type=Path)
    p.add_argument("--only", default=",".join(TREES))
    p.add_argument("--rebuild", action="store_true")
    a = p.parse_args()
    names = [n for n in a.only.split(",") if n]
    # one build at a time (machine is shared with the user)
    with ThreadPoolExecutor(1) as ex:
        for name, rc in ex.map(lambda n: build(n, a.root, a.rebuild), names):
            print(f"{name}: {'OK' if rc == 0 else f'FAILED rc={rc}'}")


if __name__ == "__main__":
    main()

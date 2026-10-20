"""Repro: Horndeski background solver keeps num_eq from a previous model (implicit SAVE of an initialised local).

fortran/eftcamb/08f_full_models/008p0_Horndeski.f90, EFTCAMBHdskSolveBackgroundEquations:
    integer :: num_eq = 5          ! initialisation in the declaration => implicit SAVE
    ...
    if (self%evolve_hubble .or. self%designerV) num_eq = 6
num_eq is never reset to 5, so after any Horndeski run with Horndeski_evolve_hubble=T (the default) or
Horndeski_designer_V=T, every later model in the same process with both flags F integrates 6 ODEs, the 6th with an
uninitialised y(6) and an undefined ydot(6).

Test: 5_hdsk_TG_1 (Horndeski_evolve_hubble = F) background in a fresh process vs in a process that first ran the same
model with Horndeski_evolve_hubble = T.  Expected: bit-identical results.

Run: /Users/marco/Software/miniforge3/bin/python hdsk-num-eq-save.py [tree]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
REPO = Path(__file__).resolve().parents[4]
INI = REPO / "fortran" / "eftcamb_test" / "parameters" / "5_hdsk_TG_1.ini"


def write_ini(tmp: Path, name: str, evolve_hubble: bool) -> Path:
    p = tmp / f"{name}.ini"
    p.write_text(
        f"DEFAULT({INI})\n"
        f"output_root = {tmp / name}\n"
        "EFTCAMB_write_background = F\n"
        "feedback_level = 0\n"
        f"Horndeski_evolve_hubble = {'T' if evolve_hubble else 'F'}\n"
    )
    return p


def num_eq_reader(camb):
    """Return a function reading the SAVEd local num_eq of the Horndeski solver from the loaded camblib.so."""
    import ctypes

    lib = Path(camb.__file__).parent / "camblib.so"
    syms = {}
    for line in subprocess.run(["nm", str(lib)], capture_output=True, text=True).stdout.splitlines():
        parts = line.split()
        if len(parts) == 3:
            syms[parts[2]] = int(parts[0], 16)
    ref = "___eftcamb_python_MOD_eftcamb_get_num_params"
    var = [k for k in syms if k.startswith("_num_eq.")]
    if not var:
        return lambda: None
    base = ctypes.cast(camb.baseconfig.camblib.__eftcamb_python_MOD_eftcamb_get_num_params, ctypes.c_void_p).value - syms[ref]
    addr = base + syms[var[0]]
    return lambda: ctypes.c_int.from_address(addr).value


def worker(tree: str, inis: list[str], out: str) -> None:
    sys.path.insert(0, tree)
    import numpy as np

    import camb

    read_num_eq = num_eq_reader(camb)
    res = []
    for ini in inis:
        cp = camb.read_ini(ini)
        r = camb.get_background(cp)
        a = np.logspace(-6, 0, 200)
        _, arr = cp.EFTCAMB.get_eft_functions(r, a)
        z = np.linspace(0, 5, 50)
        res.append(
            {
                "ini": ini,
                "H": [float(x) for x in r.hubble_parameter(z)],
                "tau0": float(r.conformal_time(0)),
                "OmegaV": [float(x) for x in arr["EFTOmegaV"]],
                "adotoa": [float(x) for x in arr["adotoa"]],
                "EFTc": [float(x) for x in arr["EFTc"]],
                "num_eq_after": read_num_eq(),
            }
        )
    Path(out).write_text(json.dumps(res))


def run(tree: str, inis: list[Path], out: Path) -> list[dict]:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    rc = subprocess.run(
        [sys.executable, __file__, "_worker", tree, json.dumps([str(i) for i in inis]), str(out)], env=env
    ).returncode
    if rc != 0:
        print(f"worker exit code {rc}")
        return []
    return json.loads(out.read_text())


def maxrel(a, b) -> float:
    import numpy as np

    a, b = np.asarray(a), np.asarray(b)
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))


def main() -> None:
    arg = sys.argv[1] if len(sys.argv) > 1 else "A0"
    tree = arg if os.path.isabs(arg) else str(ROOT / "trees" / arg)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        ini_f = write_ini(tmp, "hdsk_evolveH_F", False)
        ini_t = write_ini(tmp, "hdsk_evolveH_T", True)
        fresh = run(tree, [ini_f], tmp / "fresh.json")
        seq = run(tree, [ini_t, ini_f], tmp / "seq.json")
    if not fresh or len(seq) < 2:
        print("FAILED to run")
        sys.exit(2)
    f, s = fresh[0], seq[1]
    print(f"tree: {tree}")
    print(f"SAVEd num_eq after the run: fresh F run = {f['num_eq_after']}, T run = {seq[0]['num_eq_after']}, "
          f"F run after T run = {s['num_eq_after']}  (an evolve_hubble=F model needs 5 equations)")
    print(f"model: 5_hdsk_TG_1 (Horndeski_evolve_hubble=F); sequence: same model with evolve_hubble=T first")
    print(f"tau0       fresh = {f['tau0']:.15e}  after-T = {s['tau0']:.15e}  rel diff = {maxrel(s['tau0'], f['tau0']):.3e}")
    for key in ("H", "adotoa", "OmegaV", "EFTc"):
        print(f"{key:10s} max rel diff (after-T vs fresh) = {maxrel(s[key], f[key]):.3e}   bit identical = "
              f"{s[key] == f[key]}")
    ok = all(s[k] == f[k] for k in ("H", "adotoa", "OmegaV", "EFTc")) and s["tau0"] == f["tau0"]
    print("outputs: expected bit identical; observed:", "bit identical" if ok else "DIFFERENT")
    leak = s["num_eq_after"] is not None and s["num_eq_after"] != f["num_eq_after"]
    print("state: expected num_eq = 5 for the F model in both processes; observed:",
          f"LEAKED ({s['num_eq_after']} after a T run)" if leak else "no leak")
    ok = ok and not leak
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_worker":
        worker(sys.argv[2], json.loads(sys.argv[3]), sys.argv[4])
    else:
        main()

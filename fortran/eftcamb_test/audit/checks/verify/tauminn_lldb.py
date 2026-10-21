"""Verifier (t3-eft-turnon-before-tauminn-lmax): run 5_ADE_1 (l_max_scalar=L, no tensors) on the frozen A0 build
under lldb, breaking at Thermo_Init / GetTauStart to print the taumin actually passed (no rebuild needed).

Usage: python tauminn_lldb.py child <L> <tree> <ini> [kmax_hMpc]   (the CAMB computation; run under lldb by the parent)
       python tauminn_lldb.py <L> [tree] [ini] [kmax_hMpc]          (parent: holds a CPU slot and drives lldb)
L=0 means no CMB (transfer only).
"""

import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))
import harness  # noqa: E402


def child(lmax: int, tree: str, ini: str, kmax: str | None) -> None:
    sys.path.insert(0, str(harness.tree_path(tree)))
    import camb

    cp = camb.read_ini(ini, no_validate=True)
    if lmax > 0:
        cp.max_l = lmax
        cp.max_eta_k = 2.5 * lmax
    else:
        cp.WantCls = False
        cp.WantScalars = False
        cp.DoLensing = False
    if kmax:
        cp.Transfer.kmax = float(kmax) * cp.H0 / 100
    cp.WantTensors = False
    cp.InitPower.r = 0.0
    print("CHILD max_l", cp.max_l, "max_eta_k", cp.max_eta_k, "kmax", cp.Transfer.kmax, flush=True)
    try:
        r = camb.get_results(cp)
        print("CHILD result ok tau0", r.tau0, flush=True)
    except Exception as e:  # noqa: BLE001
        print("CHILD error", str(e).replace("\n", " | "), flush=True)


def parent(lmax: int, tree: str, ini: str, kmax: str | None) -> None:
    py = sys.executable
    cmds = [
        "settings set target.process.stop-on-exec false",
        "breakpoint set -n __results_MOD_thermo_init",
        f"command script import {HERE.parent / 'tauminn_lldb_cb.py'}",
        "breakpoint command add -F tauminn_lldb_cb.thermo_cb 1",
        "breakpoint set -n __cambmain_MOD_gettaustart",
        "breakpoint command add -F tauminn_lldb_cb.taustart_cb 2",
        "process launch -- " + f"{HERE} child {lmax} {tree} {ini} {kmax or ''}",
        "quit",
    ]
    args = ["lldb", "--batch"]
    for c in cmds:
        args += ["-o", c]
    args += ["--", py]
    env = dict(os.environ, OMP_NUM_THREADS="1")
    with harness.slot(1):
        p = subprocess.run(args, capture_output=True, text=True, env=env, timeout=600)
    lines = [ln for ln in (p.stdout + p.stderr).splitlines() if "AUDIT" in ln or "CHILD" in ln]
    # print taumin passed to Thermo_Init and the first few q passed to GetTauStart
    print("\n".join(lines[:40]))
    print("... total lines", len(lines))


if __name__ == "__main__":
    if sys.argv[1] == "child":
        child(int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5] if len(sys.argv) > 5 else None)
    else:
        a = sys.argv[1:] + [None] * 4
        parent(int(a[0]), a[1] or "A0", a[2] or str(harness.PARAMS_DIR / "5_ADE_1.ini"), a[3])

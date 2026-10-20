"""Repro (t3-gr): EFTCAMB.get_scale_evolution / get_evolution after a matter-power-only run (get_scalar_cls = F)
kills the Python process with SIGSEGV/SIGBUS.

camb/eftcamb.py get_evolution -> 11_EFTCAMB_python.f90 EFTCAMB_GetEvolution: each k mode is set up with
GetNumEqns(Ev), which only fills the scalar hierarchy sizes (lmaxg, lmaxnr, MaxlNeeded, nvar, ...) when
CP%WantScalars is true. After a transfer-only run the results carry WantScalars = F, so the OMP-private Ev keeps
stack garbage: initial() then runs `EV%Kf(1:EV%MaxlNeeded) = 1` with MaxlNeeded ~ 2^31 and writes past the stack
(lldb: EXC_BAD_ACCESS in __gaugeinterface_MOD_initial, called from geteftoutputevolutionfork). No error is raised.
(Upstream CAMB_TimeEvolution has the same unguarded GetNumEqns call; it happens not to crash in B0.)
Candidate fix (tree work/t3-gr/R8): set this%CP%WantScalars = .true. around the k loop of EFTCAMB_GetEvolution
(restore it afterwards), or raise a CAMBError when WantScalars is false.

Model: 3_RPH_MK_const_1.ini, get_scalar_cls=F, get_transfer=T (source windows removed: transfer-only + windows is a
separate upstream crash, t10a-upstream-get-background-windows-segfault.py). One subprocess per case.
Expected: return code 0 and the evolution printed. Observed: printed.
Run: /Users/marco/Software/miniforge3/bin/python t3-get-evolution-transfer-only-segfault.py [FIXED_TREE]
"""

import signal
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

INI = harness.PARAMS_DIR / "3_RPH_MK_const_1.ini"
OUT = harness.AUDIT_ROOT / "work" / "t3-gr" / "repro_getevol"

WORKER = r"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, sys.argv[1])
import camb, numpy as np
cp = camb.read_ini(sys.argv[2], no_validate=True)
want_cls = sys.argv[3] == "cls"
cp.WantCls = want_cls
cp.WantScalars = want_cls
cp.WantTensors = False
cp.InitPower.r = 0.0
cp.DoLensing = False
cp.max_l = 300
cp.WantTransfer = True
cp.SourceWindows = []
r = camb.get_results(cp)
print("get_results ok", flush=True)
names, arr = cp.EFTCAMB.get_scale_evolution(r, [0.05], np.array([0.1, 0.5, 1.0]))
print("get_scale_evolution ok: pi(k=0.05, a=0.1,0.5,1) =", np.asarray(arr["pi"])[0].tolist(), flush=True)
"""


def case(tree, mode):
    OUT.mkdir(parents=True, exist_ok=True)
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-I", "-c", WORKER, str(harness.tree_path(tree)), str(INI), mode],
                           capture_output=True, text=True, timeout=600, cwd=OUT)
    rc = p.returncode
    sig = f" ({signal.Signals(-rc).name})" if rc < 0 else ""
    last = [ln for ln in p.stdout.splitlines() if " ok" in ln][-1:]
    return f"rc={rc}{sig} {last}"


def main() -> None:
    trees = [("A0", "A0")] + ([("fixed", sys.argv[1])] if len(sys.argv) > 1 else [])
    for name, tree in trees:
        print(f"{name:6s} run with C_l (WantScalars=T)  : {case(tree, 'cls')}")
        print(f"{name:6s} transfer-only (WantScalars=F) : {case(tree, 'transfer')}")
    print("expected: rc=0 and the evolution printed in both modes; observed on A0: transfer-only killed by a signal")


if __name__ == "__main__":
    main()

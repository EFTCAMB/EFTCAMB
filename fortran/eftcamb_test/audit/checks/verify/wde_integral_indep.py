"""Independent check of finding wde-integral-unsupported-parametrizations.

Different code path from the repro: ShiftSym alpha_B (EFTflag=2, AltParEFTmodel=4, 007p7) and designer f(R),
driven directly by camb.read_ini + camb.get_results in a child process (no harness.run_model), recording the child's
exit status and whether ERROR STOP text appears, i.e. whether the Python interpreter is killed instead of a
CAMBError being raised.
"""

import os
import subprocess
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "wde_integral"
BASE = harness.PARAMS_DIR / "base_params.ini"
TREE = harness.tree_path("A0")

SHIFTSYM = f"""EFTflag = 2
AltParEFTmodel = 4
Shift_Symmetric_alphaB0 = 0.1
Shift_Symmetric_m = 2.0
Shift_Symmetric_alphaK0 = 1.0
EFT_ghost_stability = T
EFT_gradient_stability = T
EFTCAMB_write_background = F
DEFAULT({BASE})
"""
FR = f"""DEFAULT({harness.PARAMS_DIR / "4_Des_F_R_wCDM_1.ini"})
EFTCAMB_write_background = F
"""

CASES = {
    "shiftsym_wDE1": (SHIFTSYM, "EFTwDE = 1\nEFTw0 = -0.95\n"),
    "shiftsym_wDE6": (SHIFTSYM, "EFTwDE = 6\nEFTw_Pade_order_N = 0\nEFTw_Pade_order_D = 0\nEFTwa0 = 0\nEFTwN0 = -0.95\n"),
    "shiftsym_wDE7": (SHIFTSYM, "EFTwDE = 7\nEFTw_a0 = -0.95\n"),
    "shiftsym_wDE11": (SHIFTSYM, "EFTwDE = 11\nEFTw0 = -0.95\nEFTw1 = 0.0\n"),
    "fR_wDE11_lam0.5": (FR, "EFTwDE = 11\nEFTw0 = -0.9\nEFTw1 = 0.5\n"),
}

CHILD = textwrap.dedent("""
    import sys
    sys.path.insert(0, sys.argv[1])
    import camb
    try:
        pars = camb.read_ini(sys.argv[2])
        r = camb.get_results(pars)
        print("RESULT ok sigma8=%.6f" % r.get_sigma8_0())
    except Exception as e:
        print("RESULT python-exception", type(e).__name__, str(e)[:120])
""")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    for name, (head, extra) in CASES.items():
        ini = OUT / f"{name}.ini"
        ini.write_text(extra + head)
        with harness.slot(1):
            p = subprocess.run([sys.executable, "-c", CHILD, str(TREE), str(ini)], capture_output=True, text=True,
                               env=env, cwd=OUT, timeout=900)
        txt = p.stdout + p.stderr
        res = [ln for ln in txt.splitlines() if ln.startswith("RESULT")]
        stop = [ln.strip() for ln in txt.splitlines() if "not implemented" in ln or "ERROR STOP" in ln]
        print(f"{name:18s} rc={p.returncode:3d} result={res[-1] if res else 'NONE (interpreter died)'} "
              f"| {' / '.join(stop[:2])}")


if __name__ == "__main__":
    main()

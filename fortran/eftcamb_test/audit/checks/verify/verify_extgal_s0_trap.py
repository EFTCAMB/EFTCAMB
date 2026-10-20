"""Adversarial verifier (t3-extgal-gr-point-background): locate the first floating-point exception at S = 0.

Runs the Extended Galileon at S = 0 (skip stability) on an AD tree copy (-ffpe-trap=invalid,zero,overflow) that only
has the F-002 nested-if fix (work/t1-crash/chk/AD, read-only), and prints the backtrace lines from the log.
"""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

TREE = harness.AUDIT_ROOT / "work/t1-crash/chk/AD"
INI = REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"
OUT = harness.AUDIT_ROOT / "work/verify/extgal_gr/trap"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "200", "get_transfer": "F",
        "num_redshiftwindows": "0", "EFTCAMB_skip_stability": "T"}
for name, ov in (("S0", {"ExtendedGalileon_S": "0", "ExtendedGalileon_q": "100"}),):
    st = harness.run_model(TREE, INI, OUT / f"{name}.npz", overrides={**FAST, **ov}, threads=1, timeout=540)
    print(name, st["status"])
    log = Path(st["log"]).read_text(errors="replace").splitlines()
    for ln in log:
        if "signal" in ln.lower() or ".f90:" in ln or "Floating" in ln:
            print("   ", ln.strip())

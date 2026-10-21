"""Repro (t4a-param): quintic-spline parametrization with the accepted option EFT*_Spline_Order = 3
(init_func_from_file accepts 3 or 5) stops the program in third_derivative/fourth_derivative
(04p019_spline5_parametrizations_1D.f90:541-544, 597-600: MpiStop -> ERROR STOP).  Pure EFT always evaluates
EFTOmegaPPP/PPPP in the background, so any Pure-EFT Omega spline of order 3 kills the Python interpreter.
The message is also garbled: format '(a,a,a,I3)' is given (string, string, integer).

Model: 2_PEFT_Omega_lin_1.ini (Omega = a) with Omega given as an order-3 spline of exactly linear nodes;
expected: runs and reproduces the reference (order 5 with the same nodes does).

usage: /Users/marco/Software/miniforge3/bin/python spline5-order3-errorstop.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t4a-param" / "repro" / Path(__file__).stem
INI = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"
X6 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    ref = WORK / "ref.npz"
    st_ref = harness.run_model("A0", INI, ref)
    for order in (5, 3):
        ov = {"PureEFTmodelOmega": "10", "EFTOmega_Spline_Pixels": "6", "EFTOmega_Spline_Order": str(order)}
        ov.update({f"EFTOmegax{i + 1}": str(x) for i, x in enumerate(X6)})
        ov.update({f"EFTOmegav{i + 1}": str(x) for i, x in enumerate(X6)})
        out = WORK / f"order{order}.npz"
        st = harness.run_model("A0", INI, out, overrides=ov)
        msg = st["status"]
        if st["status"] == "ok" and st_ref["status"] == "ok":
            msg += f", worst Cl rel diff vs Omega=a reference = {harness.compare(out, ref)['worst_cl']:.2e}"
        else:
            tail = [repr(t) for t in (st.get("log_tail") or []) if t.strip()]
            msg += " | " + " / ".join(tail[:1] + tail[-3:])
        print(f"Spline_Order = {order}: expected ok (~0 difference) | observed: {msg}")


if __name__ == "__main__":
    main()

"""Repro (t4a-param): designer f(R) (EFTflag=3, DesignerEFTmodel=1) with the documented EFTwDE = 6..11 options.

Every case below describes the SAME constant w = -0.9 as the test model 4_Des_F_R_wCDM_1.ini (EFTwDE = 1), so the
expected result is identical spectra.  Designer f(R) needs integral() (grhov_t = grhov*integral(a)):
  * Pade (6), Fourier (7), exp2 (11): integral() calls MpiStop -> ERROR STOP kills the Python interpreter
  * cubic spline (9), quintic spline (10): integral() has the wrong convention / NaN -> background solver fails
  * steplog (8) with v1 = v2: works (control)

usage: /Users/marco/Software/miniforge3/bin/python wde-integral-unsupported-parametrizations.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t4a-param" / "repro" / Path(__file__).stem
INI = harness.PARAMS_DIR / "4_Des_F_R_wCDM_1.ini"
X6 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
LOGX9 = [10.0 ** e for e in range(-8, 1)]

CASES = {
    "EFTwDE=1 constant (reference)": {},
    "EFTwDE=6 Pade N=D=0": {"EFTwDE": 6, "EFTw_Pade_order_N": 0, "EFTw_Pade_order_D": 0, "EFTwa0": 0, "EFTwN0": -0.9},
    "EFTwDE=7 Fourier a0 only": {"EFTwDE": 7, "EFTw_a0": -0.9},
    "EFTwDE=8 steplog v1=v2": {"EFTwDE": 8, "EFTw_v1": -0.9, "EFTw_v2": -0.9, "EFTw_at": 0.5, "EFTw_delta": 0.3},
    "EFTwDE=9 cubic spline": {"EFTwDE": 9, "EFTw_Spline_Pixels": 6, "EFTw_null_value": -0.9,
                              **{f"EFTwx{i + 1}": x for i, x in enumerate(X6)},
                              **{f"EFTwv{i + 1}": -0.9 for i in range(6)}},
    "EFTwDE=10 quintic spline logx": {"EFTwDE": 10, "EFTw_Spline_Pixels": 9, "EFTw_null_value": -0.9,
                                      "EFTw_Logx_Flag": "T", **{f"EFTwx{i + 1}": x for i, x in enumerate(LOGX9)},
                                      **{f"EFTwv{i + 1}": -0.9 for i in range(9)}},
    "EFTwDE=11 exp2 lambda=0": {"EFTwDE": 11, "EFTw0": -0.9, "EFTw1": 0.0},
}


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    ref = None
    for i, (name, ov) in enumerate(CASES.items()):
        out = WORK / f"case{i}.npz"
        st = harness.run_model("A0", INI, out, overrides={k: str(v) for k, v in ov.items()}, threads=2)
        tail = [line for line in (st.get("log_tail") or []) if line.strip()][-3:]
        msg = st["status"]
        if st["status"] == "ok":
            if ref is None:
                ref = out
                msg += " (reference)"
            else:
                c = harness.compare(out, ref)
                msg += f", worst Cl rel diff vs reference = {c['worst_cl']:.2e}"
        else:
            killed = any("ERROR STOP" in line for line in (st.get("log_tail") or []))
            msg += " | interpreter killed by ERROR STOP" if killed else " | CAMBError raised"
            msg += " | " + " / ".join(t.strip() for t in tail)
        print(f"{name:32s} expected: ok, identical to reference | observed: {msg}")


if __name__ == "__main__":
    main()

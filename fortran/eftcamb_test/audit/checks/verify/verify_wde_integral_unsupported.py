"""Verifier: ShiftSym alphaB (007p7) with EFTwDE = 6/7/11, and designer f(R) with Pade N=1 (non-trivial)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "verify" / "wde-integral-unsupported"
SS = harness.AUDIT_ROOT / "work" / "verify" / "wde-integral-unsupported" / "shiftsym.ini"
FR = harness.PARAMS_DIR / "4_Des_F_R_wCDM_1.ini"
CASES = [
    ("SS EFTwDE=1 ref", SS, {}),
    ("SS EFTwDE=7 Fourier", SS, {"EFTwDE": 7, "EFTw_a0": -0.95}),
    ("SS EFTwDE=11 exp2", SS, {"EFTwDE": 11, "EFTw0": -0.95, "EFTw1": 0.0}),
    ("SS EFTwDE=6 Pade", SS, {"EFTwDE": 6, "EFTw_Pade_order_N": 0, "EFTw_Pade_order_D": 0, "EFTwa0": 0,
                              "EFTwN0": -0.95}),
    ("fR EFTwDE=6 Pade N=1x", FR, {"EFTwDE": 6, "EFTw_Pade_order_N": 1, "EFTw_Pade_order_D": 0, "EFTwa0": 1,
                                  "EFTwN0": -0.9, "EFTwN1": 0.05}),
]
WORK.mkdir(parents=True, exist_ok=True)
for i, (name, ini, ov) in enumerate(CASES):
    st = harness.run_model("A0", ini, WORK / f"c{i}.npz", overrides={k: str(v) for k, v in ov.items()}, threads=1)
    tail = [t.strip() for t in (st.get("log_tail") or []) if t.strip()][-3:]
    print(f"{name:24s} {st['status']} | {' / '.join(tail) if st['status'] != 'ok' else ''}")

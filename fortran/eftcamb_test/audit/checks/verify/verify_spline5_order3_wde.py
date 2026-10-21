"""Verifier (spline5-order3-errorstop): independent path. Use the quintic-spline type for Pure EFT w_DE
(EFTwDE=10, constant nodes w=-1) instead of Omega; Pure EFT evaluates PureEFTwDE%third_derivative in the
background. Order 5 should reproduce EFTwDE=0; order 3 is accepted by init but is claimed to ERROR STOP.
Also run on B1 (pre-rebase) to see whether the behaviour predates the rebase."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "verify" / "spline5-order3"
INI = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"
X = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    ref = WORK / "ref_wLCDM.npz"
    st = harness.run_model("A0", INI, ref, threads=1)
    print("ref", st["status"])
    for tree in ("A0", "B1"):
        for order in (5, 3):
            ov = {"EFTwDE": "10", "EFTw_Spline_Pixels": "6", "EFTw_Spline_Order": str(order)}
            ov.update({f"EFTwx{i + 1}": str(x) for i, x in enumerate(X)})
            ov.update({f"EFTwv{i + 1}": "-1.0" for i in range(6)})
            out = WORK / f"{tree}_wspline_o{order}.npz"
            s = harness.run_model(tree, INI, out, overrides=ov, threads=1)
            msg = s["status"]
            if msg == "ok":
                m = harness.load(out)["meta"]
                msg += f" model={m.get('eft_model_name') if isinstance(m, dict) else ''}"
                if tree == "A0":
                    msg += f" worst_cl vs wDE=0: {harness.compare(out, ref)['worst_cl']:.2e}"
            else:
                tail = [t for t in (s.get("log_tail") or []) if t.strip()]
                msg += " | " + " / ".join(repr(t) for t in tail[-4:])
            print(f"{tree} order={order}: {msg}")


if __name__ == "__main__":
    main()

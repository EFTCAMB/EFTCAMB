"""Verifier (spline5-order3-errorstop): CAMB-level attribution with other call sites than the reporter's.
  g1_o3   : Pure EFT, Gamma1 = spline5 order 3 (only value/d1/d2 used)            -> expect ok (order 3 usable here)
  g2_o3   : Pure EFT, Gamma2 = spline5 order 3 (Gamma2PPP at 007p1:1556)           -> expect ERROR STOP if claim holds
  g2_o5   : same nodes, order 5                                                    -> expect ok
  om_o3_B1: reporter's Omega order-3 config on the pre-rebase tree B1              -> inherited or rebase regression?
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "verify" / "spline5_o3"
INI = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"
X6 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]


def spline(name, order, vals):
    ov = {f"{name}_Spline_Pixels": "6", f"{name}_Spline_Order": str(order)}
    ov.update({f"{name}x{i + 1}": str(x) for i, x in enumerate(X6)})
    ov.update({f"{name}v{i + 1}": str(v) for i, v in enumerate(vals)})
    return ov


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    small = [0.01 * x**2 for x in X6]
    cases = {
        "g1_o3": ("A0", {"PureEFTmodelGamma1": "10", **spline("EFTGamma1", 3, small)}),
        "g2_o3": ("A0", {"PureEFTmodelGamma2": "10", **spline("EFTGamma2", 3, small)}),
        "g2_o5": ("A0", {"PureEFTmodelGamma2": "10", **spline("EFTGamma2", 5, small)}),
        "om_o3_B1": ("B1", {"PureEFTmodelOmega": "10", **spline("EFTOmega", 3, X6)}),
        "om_o5_B1": ("B1", {"PureEFTmodelOmega": "10", **spline("EFTOmega", 5, X6)}),
    }
    ref = WORK / "A0_ref.npz"
    if not ref.exists():
        harness.run_model("A0", INI, ref, threads=1)
    for name, (tree, ov) in cases.items():
        out = WORK / f"{name}.npz"
        st = harness.run_model(tree, INI, out, overrides=ov, threads=1)
        msg = st["status"]
        if st["status"] == "ok":
            r = WORK / f"{tree}_ref.npz"
            if not r.exists():
                harness.run_model(tree, INI, r, threads=1)
            msg += f", worst Cl rel diff vs {tree} Omega=a reference = {harness.compare(out, r)['worst_cl']:.2e}"
        else:
            tail = [repr(t) for t in (st.get("log_tail") or []) if t.strip()]
            msg += " | " + " / ".join(tail[-5:])
        print(f"{name:9s} [{tree}]: {msg}", flush=True)


if __name__ == "__main__":
    main()

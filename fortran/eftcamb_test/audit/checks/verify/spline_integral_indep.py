"""Verifier: designer f(R) (4_Des_F_R_wCDM_1) with constant w=-0.9 given as spline (EFTwDE=9/10) vs EFTwDE=1, full
CAMB run through the Python wrapper on A0; plus numpy re-implementation of the cubic yint recursion."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = Path(harness.AUDIT_ROOT) / "work" / "verify" / "spline-integral"
INI = harness.PARAMS_DIR / "4_Des_F_R_wCDM_1.ini"
w = -0.9


def nodes(xs, flag, extra=None):
    d = {"EFTwDE": str(flag), "EFTw_Spline_Pixels": str(len(xs)), "EFTw_null_value": str(w)}
    for i, x in enumerate(xs):
        d[f"EFTwx{i + 1}"] = repr(x)
        d[f"EFTwv{i + 1}"] = repr(w)
    d.update(extra or {})
    return d


cases = {
    "ref_w1": {"EFTwDE": "1", "EFTw0": str(w)},
    "cubic_x0": nodes([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], 9),
    "cubic_x1e-8": nodes([1e-8, 0.2, 0.4, 0.6, 0.8, 1.0], 9),
    "cubic_x1e-8_to_1.1": nodes([1e-8, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1], 9),
    "quintic_log": nodes([10.0 ** e for e in range(-8, 1)], 10, {"EFTw_Logx_Flag": "T"}),
}
if __name__ == "__main__":
    res = {}
    for name, ov in cases.items():
        st = harness.run_model("A0", INI, OUT / f"{name}.npz", overrides=ov, threads=1, timeout=900)
        res[name] = {"status": st["status"], "tail": " | ".join(map(str, st.get("log_tail") or []))[-400:]}
        print(name, st["status"], " | ".join(map(str, st.get("log_tail") or []))[-400:])
    ref = OUT / "ref_w1.npz"
    for name in cases:
        p = OUT / f"{name}.npz"
        if name != "ref_w1" and p.exists() and ref.exists():
            c = harness.compare(ref, p)
            print(name, "worst_cl", c["worst_cl"])
    # numpy re-implementation of LinearIntepolateFunction1DInitIntegration for constant w on nodes 1e-8,0.2..1
    x = np.array([1e-8, 0.2, 0.4, 0.6, 0.8, 1.0])
    y = np.full_like(x, w)
    yint = np.ones_like(x)
    for i in range(len(x) - 2, -1, -1):
        b = (x[i] * y[i + 1] - x[i + 1] * y[i]) / (x[i] - x[i + 1])
        yint[i] = yint[i + 1] * (x[i] / x[i + 1]) ** (-3 * (1 + b)) * np.exp(-3 * (y[i] - y[i + 1]))
    print("numpy yint nodes", yint)
    print("exact a^1.7 at nodes", x ** 1.7, " a^-0.3 (no a^2 factor)", x ** -0.3)
    json.dump(res, open(OUT / "summary.json", "w"), indent=1)

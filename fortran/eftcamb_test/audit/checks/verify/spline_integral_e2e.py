"""Verifier (spline-integral-convention): end-to-end designer f(R) with constant w=-0.9 given as EFTwDE=1 vs 9/10."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "spline_integral_e2e"
OUT.mkdir(parents=True, exist_ok=True)
INI = harness.model_inis("4_Des_F_R_wCDM_1.ini")[0]
w = -0.9


def spline(flag, xs, extra=None):
    ov = {"EFTwDE": str(flag), "EFTw_Spline_Pixels": str(len(xs)), "EFTw_null_value": str(w)}
    for i, x in enumerate(xs):
        ov[f"EFTwx{i + 1}"] = repr(x)
        ov[f"EFTwv{i + 1}"] = str(w)
    ov.update(extra or {})
    return ov


cases = {
    "const": {"EFTwDE": "1", "EFTw0": str(w)},
    "cubic_lin": spline(9, [1e-8, 0.2, 0.4, 0.6, 0.8, 1.0]),
    "cubic_past1": spline(9, [1e-8, 0.25, 0.5, 0.75, 1.0, 1.25]),
    "quintic_log": spline(10, [10.0**e for e in range(-8, 1)], {"EFTw_Logx_Flag": "T"}),
    "quintic_log_past1": spline(10, [10.0**e for e in range(-8, 1)] + [1.5], {"EFTw_Logx_Flag": "T"}),
}
res = {}
for name, ov in cases.items():
    st = harness.run_model("A0", INI, OUT / f"{name}.npz", overrides=ov, threads=1, timeout=600)
    res[name] = {"status": st["status"], "tail": " | ".join(map(str, st.get("log_tail") or []))[-400:]}
    print(name, st["status"], " | ".join(map(str, st.get("log_tail") or []))[-400:])
ref = OUT / "const.npz"
if ref.exists():
    r = harness.load(ref)
    a = r["eft_a_grid"]
    for name in cases:
        p = OUT / f"{name}.npz"
        if name == "const" or not p.exists():
            continue
        d = harness.load(p)
        c = harness.compare(r, d)
        print(name, "worst_cl", c["worst_cl"], "worst_pk", c["worst_pk"])
        for fld in ("eft_grhov_t", "eft_EFTOmegaV"):
            if fld in d:
                m = a > 1e-4
                print("  ", fld, "max rel diff a>1e-4:", np.nanmax(np.abs(d[fld][m] / r[fld][m] - 1)))
json.dump(res, open(OUT / "status.json", "w"), indent=1)

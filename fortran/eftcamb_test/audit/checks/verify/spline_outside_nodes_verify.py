"""Independent check for spline-derivatives-outside-nodes: same constant w=-1 spline, but nodes covering the whole
a range queried (no 'outside' region) vs nodes 0.01..0.9; and null_value=0 cases. Compare each to EFTwDE=0."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness

W = harness.AUDIT_ROOT / "work" / "verify" / "spline-derivatives-outside-nodes"
W.mkdir(parents=True, exist_ok=True)
ini = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"


def keys(nodes, null):
    d = {"EFTwDE": "9", "EFTw_Spline_Pixels": str(len(nodes)), "EFTw_null_value": str(null)}
    d.update({f"EFTwx{i + 1}": repr(x) for i, x in enumerate(nodes)})
    d.update({f"EFTwv{i + 1}": "-1.0" for i in range(len(nodes))})
    return d


ref = W / "ref.npz"
if not ref.exists():
    print(harness.run_model("A0", ini, ref, threads=1)["status"])
cases = {
    "narrow_null-1": keys([0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9], -1.0),
    "wide_null-1": keys([1e-12, 1e-6, 1e-3, 0.01, 0.1, 0.3, 0.5, 0.9, 1.5, 2.0], -1.0),
}
for name, ov in cases.items():
    out = W / f"{name}.npz"
    st = harness.run_model("A0", ini, out, overrides=ov, threads=1)
    print(name, st["status"])
    if st["status"] != "ok":
        print(st.get("log_tail"))
        continue
    a, r = harness.load(out), harness.load(ref)
    c = harness.compare(a, r)
    print(f"  worst_cl {c['worst_cl']:.3e} worst_pk {c['worst_pk']:.3e} eft top:",
          sorted(c["eft"].items(), key=lambda kv: -kv[1])[:3])
    ag = a["eft_a_grid"]
    for f in ("eft_EFTLambdadot",):
        if f in a:
            d = np.abs(a[f] - r[f])
            i = np.argmax(d)
            print(f"  {f}: max|diff| {d[i]:.3e} at a={ag[i]:.3e}; diff at a<0.01 max {d[ag < 0.01].max():.3e},"
                  f" a>0.9 max {d[ag > 0.9].max():.3e}, 0.02<a<0.85 max {d[(ag > 0.02) & (ag < 0.85)].max():.3e}")

"""Independent check for spline-derivatives-outside-nodes: if the null_value derivatives are the cause, the same w=-1
spline with nodes covering the whole integration range (so no point is ever outside) must reproduce EFTwDE=0, while
the narrow-node version must not. Also inspect EFTLambdadot(a) where it differs."""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "verify" / "spline_outside"
W.mkdir(parents=True, exist_ok=True)
ini = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"


def keys(flag, nodes):
    d = {"EFTwDE": flag, "EFTw_Spline_Pixels": len(nodes), "EFTw_null_value": -1.0}
    if flag == 10:
        d.update({"EFTw_Spline_Order": 5, "EFTw_Logx_Flag": "F", "EFTw_Logy_Flag": "F"})
    d.update({f"EFTwx{i + 1}": x for i, x in enumerate(nodes)})
    d.update({f"EFTwv{i + 1}": -1.0 for i in range(len(nodes))})
    return {k: str(v) for k, v in d.items()}


wide = [1e-12, 1e-6, 1e-3, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.5]
narrow = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
tree = sys.argv[1] if len(sys.argv) > 1 else "A0"
runs = {"ref": {}}
for flag in (9, 10):
    runs[f"wide{flag}"] = keys(flag, wide)
    runs[f"narrow{flag}"] = keys(flag, narrow)
# null_value = 0 variant: value outside is 0 (w=0 outside!) - not used. Instead narrow nodes with null_value -1 is the bug.
res = {}
for name, ov in runs.items():
    out = W / f"{tree}_{name}.npz"
    st = harness.run_model(tree, ini, out, overrides=ov, threads=1)
    print(name, st["status"])
    res[name] = harness.load(out)
ref = res["ref"]
for name in runs:
    if name == "ref":
        continue
    c = harness.compare(res[name], ref)
    print(f"{tree} {name}: worst_cl {c['worst_cl']:.3e} worst_pk {c['worst_pk']:.3e} sigma8 {c['sigma8']:.3e} "
          f"bit_identical {c['bit_identical']} EFTLambdadot {c['eft'].get('EFTLambdadot', float('nan')):.3e}")
a = ref["eft_a_grid"]
d = res["narrow10"]["eft_EFTLambdadot"] - ref["eft_EFTLambdadot"]
for ai in (1e-3, 5e-3, 0.02, 0.5, 0.95, 1.0):
    i = np.argmin(abs(a - ai))
    print(f"a={a[i]:.4g} Lambdadot ref {ref['eft_EFTLambdadot'][i]:.6e} narrow10 {res['narrow10']['eft_EFTLambdadot'][i]:.6e}"
          f" diff {d[i]:.3e}")

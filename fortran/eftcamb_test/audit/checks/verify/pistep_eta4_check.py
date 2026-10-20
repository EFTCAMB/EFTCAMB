"""Verifier (pi-rk45-step-instability): secondary models with a constant step cap and B1."""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "pistep_out"
CAP = harness.AUDIT_ROOT / "work" / "verify" / "pistep_cap"
B2 = {"accuracy_boost": "2", "l_accuracy_boost": "2", "l_sample_boost": "2"}
NW = {"num_redshiftwindows": "0"}
C2 = {"AUDIT_VFY_MAXSTEP": sys.argv[1] if len(sys.argv) > 1 else "2"}


def job(a):
    tag, tree, model, over, env = a
    out = OUT / f"{model}_{tag}.npz"
    if not out.exists():
        st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=over, threads=1,
                               extra_env=env, timeout=1200)
        if st["status"] != "ok":
            raise SystemExit(f"{model} {tag}: {st['status']}")
    return (model, tag), harness.load(out)["cl_unlensed_scalar"][:, 0]


jobs = []
for m in ("5_Horava_eta_4", "3_RPH_MB_OmegaDE_1"):
    jobs += [("A0def", "A0", m, NW, {}), ("A0b2", "A0", m, {**NW, **B2}, {}), (f"cap{C2['AUDIT_VFY_MAXSTEP']}", CAP, m, NW, C2)]
    if C2["AUDIT_VFY_MAXSTEP"] == "2":
        jobs += [("cap2b2", CAP, m, {**NW, **B2}, C2), ("B1def", "B1", m, NW, {}), ("B1b2", "B1", m, {**NW, **B2}, {})]
with ThreadPoolExecutor(4) as ex:
    cl = dict(ex.map(job, jobs))
ells = np.arange(1500, 2501)
for m in ("5_Horava_eta_4", "3_RPH_MB_OmegaDE_1"):
    cx = f"cap{C2['AUDIT_VFY_MAXSTEP']}"
    pairs = (("A0def", "A0b2"), (cx, "A0b2"))
    if cx == "cap2":
        pairs += (("cap2", "cap2b2"), ("cap2b2", "A0b2"), ("B1def", "B1b2"))
    for x, y in pairs:
        r = cl[(m, x)][ells] / cl[(m, y)][ells] - 1
        i = int(np.argmax(np.abs(r)))
        print(f"{m:20s} {x:7s}/{y:7s} max {r[i]:+.3e} at ell {ells[i]}")

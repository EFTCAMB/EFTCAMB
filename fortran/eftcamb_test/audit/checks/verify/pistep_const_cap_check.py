"""Verifier (pi-rk45-step-instability): independent check with a constant RK45 step cap and with B1 (dverk)."""
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "pistep_out"
REF = harness.AUDIT_ROOT / "work" / "t10b-convergence" / "repro_pi_rk45"
CAP = harness.AUDIT_ROOT / "work" / "verify" / "pistep_cap"
B2 = {"accuracy_boost": "2", "l_accuracy_boost": "2", "l_sample_boost": "2"}
NOWIN = {"num_redshiftwindows": "0"}
M = "5_Horava_SolSyst_2"


def job(a):
    tag, tree, model, over, env = a
    out = OUT / f"{model}_{tag}.npz"
    t0 = time.time()
    if not out.exists():
        st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=over, threads=1,
                               extra_env=env, timeout=1200)
        if st["status"] != "ok":
            raise SystemExit(f"{model} {tag}: {st['status']}")
    return (model, tag), (harness.load(out)["cl_unlensed_scalar"][:, 0], time.time() - t0)


jobs = [("cap_none", CAP, M, NOWIN, {}),
        ("cap2", CAP, M, NOWIN, {"AUDIT_VFY_MAXSTEP": "2"}),
        ("cap2_win", CAP, M, {}, {"AUDIT_VFY_MAXSTEP": "2"}),
        ("cap0.5", CAP, M, NOWIN, {"AUDIT_VFY_MAXSTEP": "0.5"}),
        ("cap5", CAP, M, NOWIN, {"AUDIT_VFY_MAXSTEP": "5"}),
        ("B1def", "B1", M, NOWIN, {}), ("B1b2", "B1", M, {**NOWIN, **B2}, {}), ("B1win", "B1", M, {}, {})]
with ThreadPoolExecutor(4) as ex:
    res = dict(ex.map(job, jobs))
cl = {k: v[0] for k, v in res.items()}
for k, v in res.items():
    print("time", k, round(v[1], 1))
ref = {t: harness.load(REF / f"{M}_{t}.npz")["cl_unlensed_scalar"][:, 0] for t in ("def", "b2", "win")}
ells = np.arange(1500, 2501)


def rep(lab, x, y):
    r = x[ells] / y[ells] - 1
    i = int(np.argmax(np.abs(r)))
    print(f"{lab:36s} max {r[i]:+.3e} at ell {ells[i]}; ell2000 {r[500]:+.3e}")


rep("cap tree, no env vs A0 def", cl[(M, "cap_none")], ref["def"])
print("bit identical:", np.array_equal(cl[(M, "cap_none")], ref["def"]))
for t in ("cap5", "cap2", "cap0.5"):
    rep(f"{t} vs A0 boost2", cl[(M, t)], ref["b2"])
rep("cap2+windows vs cap2", cl[(M, "cap2_win")], cl[(M, "cap2")])
rep("cap0.5 vs cap2", cl[(M, "cap0.5")], cl[(M, "cap2")])
rep("B1 def vs B1 boost2", cl[(M, "B1def")], cl[(M, "B1b2")])
rep("B1 def+windows vs B1 def", cl[(M, "B1win")], cl[(M, "B1def")])
rep("B1 boost2 vs A0 boost2", cl[(M, "B1b2")], ref["b2"])

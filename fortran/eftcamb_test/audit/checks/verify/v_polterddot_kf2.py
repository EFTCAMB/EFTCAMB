"""Verifier for t3-polterddot-curvature-gr: A0 vs B0 vs R3 (A0 + maxz50 + Kf(2) fix), lensed Cls, omk=+-0.1."""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "polterddot_kf2"
R3 = harness.AUDIT_ROOT / "work" / "t3-gr" / "R3"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
TREES = {"A0": "A0", "B0": "B0", "R3": str(R3)}


def job(a):
    name, omk = a
    ov = {"omk": omk, "num_redshiftwindows": "0", "get_tensor_cls": "F"}
    return harness.run_model(TREES[name], GR, OUT / f"{name}_omk{omk}.npz", overrides=ov,
                             python_set={"Reion.max_redshift": 50.0}, threads=1)


omks = ("-0.1", "0.1")
with ThreadPoolExecutor(3) as ex:
    for st in ex.map(job, [(t, o) for o in omks for t in TREES]):
        if st["status"] != "ok":
            print(st)
for o in omks:
    d = {t: harness.load(OUT / f"{t}_omk{o}.npz") for t in TREES}
    for key in ("cl_unlensed_scalar", "cl_total"):
        b = d["B0"][key]
        for t in ("A0", "R3"):
            a = d[t][key]
            r = np.abs(a[2:2000, 0] / b[2:2000, 0] - 1)
            ee = np.abs(a[2:2000, 1] / b[2:2000, 1] - 1)
            print(f"omk={o:>5s} {key:20s} {t} vs B0: max|dTT/TT| l<30 {r[:28].max():.3e} (l={2+np.argmax(r[:28])})"
                  f"  l>=100 {r[98:].max():.3e}  max|dEE/EE| {ee.max():.3e}")

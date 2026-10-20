"""Compare A0, P (pigddot Kf fix), E (E2ddot Kf fix) against B0 for GR at several omk, max_redshift=50."""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "verify" / "polter"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
TREES = {"A0": "A0", "B0": "B0", "P": W / "P", "E": W / "E"}
OMKS = ("-0.05", "0.05", "0")


def job(a):
    name, omk = a
    ov = {"omk": omk, "num_redshiftwindows": "0", "get_tensor_cls": "F", "l_max_scalar": "1500", "do_lensing": "F"}
    out = W / "runs" / f"{name}_omk{omk}.npz"
    if out.exists():
        return {"status": "ok"}
    return harness.run_model(TREES[name], GR, out, overrides=ov, python_set={"Reion.max_redshift": 50.0}, threads=1)


with ThreadPoolExecutor(4) as ex:
    for st in ex.map(job, [(n, o) for o in OMKS for n in TREES]):
        if st["status"] != "ok":
            print(st)
for o in OMKS:
    b = harness.load(W / "runs" / f"B0_omk{o}.npz")["cl_unlensed_scalar"]
    for n in ("A0", "P", "E"):
        a = harness.load(W / "runs" / f"{n}_omk{o}.npz")["cl_unlensed_scalar"]
        r = np.abs(a[2:1500, 0] / b[2:1500, 0] - 1)
        e = np.abs(a[2:1500, 1] / b[2:1500, 1] - 1)
        te = np.abs(a[2:1500, 3] - b[2:1500, 3]) / np.sqrt(b[2:1500, 0] * b[2:1500, 1])
        print(f"omk={o:>6s} {n}: TT l<30 {r[:28].max():.3e} (l={2 + r[:28].argmax()})  TT l>=100 {r[98:].max():.3e}"
              f"  EE max {e.max():.3e} (l={2 + e.argmax()})  TE {te.max():.3e}")

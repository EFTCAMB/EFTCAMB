"""Verifier for t3-eft-turnon-before-tauminn-lmax: scan turn-on/lmax, B1 behaviour, and sensitivity of Cls to turn-on."""
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

INI = harness.PARAMS_DIR / "5_ADE_1.ini"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "tauminn_lmax"
TREE = sys.argv[1] if len(sys.argv) > 1 else "A0"
CASES = [("B1", 2000, "1.e-10"), ("A0", 2000, "1.1e-10"), ("A0", 2500, "1.e-10"), ("A0", 2500, "1.1e-10"),
         ("A0", 1500, "1.e-10"), ("A0", 3000, "1.e-10")]
if TREE != "A0":
    CASES = [(TREE, 2000, "1.e-10"), (TREE, 2000, "1.1e-10"), (TREE, 2500, "1.e-10")]


def job(c):
    tree, lmax, ton = c
    tag = Path(tree).name
    out = OUT / f"{tag}_l{lmax}_t{ton}.npz"
    st = harness.run_model(tree, INI, out, overrides={"l_max_scalar": str(lmax), "get_tensor_cls": "F",
                                                      "EFTCAMB_turn_on_time": ton}, threads=1, timeout=900)
    msg = re.findall(r"(EFTCAMB starts before[^\n]*)", out.with_suffix(".log").read_text(errors="replace"))
    return c, st["status"], msg[:1], out


res = {}
with ThreadPoolExecutor(3) as ex:
    for c, s, m, out in ex.map(job, CASES):
        print(c, s, m)
        res[c] = out if s == "ok" else None


def cmp(a, b, lmax=1800):
    if res.get(a) is None or res.get(b) is None:
        return
    A, B = harness.load(res[a]), harness.load(res[b])
    for name in ("cl_total",):
        n = min(lmax, len(A[name]), len(B[name])); x, y = A[name][2:n], B[name][2:n]
        print(a, "vs", b, name, "max rel TT", np.max(np.abs(x[:, 0] / y[:, 0] - 1)), "EE", np.max(np.abs(x[:, 1] / y[:, 1] - 1)))


t0 = CASES[0][0]
cmp(("A0", 2500, "1.e-10"), ("A0", 2500, "1.1e-10"))
cmp(("A0", 2000, "1.1e-10"), ("A0", 2500, "1.1e-10"))
if TREE != "A0":
    cmp((TREE, 2000, "1.e-10"), (TREE, 2000, "1.1e-10"))
    cmp((TREE, 2500, "1.e-10"), ("A0", 2500, "1.e-10")) if False else None

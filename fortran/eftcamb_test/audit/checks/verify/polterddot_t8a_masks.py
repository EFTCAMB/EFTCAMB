"""Verifier (t8a polterddot-nonflat-gr): which Kf(2) factors are needed to reproduce CAMB.

Reuses the diagnostic tree built by polterK_attrib.py (env AUDIT_POLTERK bit1 = Kf(2) on octgdot in pigddot, bit2 =
Kf(2) on Edot(3) in E2ddot). Bit1 alone is the "hierarchy-consistent" form (denlk2(2)); bits 1+2 equal CAMB's
original polterddot. Same settings and B0 reference as polterddot_t8a_fixtree.py (runs/B0_GR_omk*.npz).
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402
from polterddot_t8a_fixtree import GR, OV, WORK, stats  # noqa: E402

DIAG = harness.AUDIT_ROOT / "work" / "verify" / "polterK" / "attrib_tree"
OMKS = ("-0.05", "0.05", "-0.1")


def job(a):
    omk, mask = a
    return harness.run_model(str(DIAG), GR, WORK / "runs" / f"MASK{mask}_GR_omk{omk}.npz", overrides={**OV, "omk": omk},
                             python_set={"Reion.max_redshift": 50.0}, threads=1,
                             extra_env={"AUDIT_POLTERK": str(mask)})


with ThreadPoolExecutor(4) as ex:
    for st in ex.map(job, [(o, m) for o in OMKS for m in (0, 1, 2, 3)]):
        if st["status"] != "ok":
            print("FAILED", st["out"], st["status"])
L = lambda n: harness.load(WORK / "runs" / f"{n}.npz")  # noqa: E731
labels = {0: "no Kf (= A0)", 1: "Kf on octg only", 2: "Kf on Edot3 only", 3: "Kf on both (CAMB)"}
for o in OMKS:
    b0 = L(f"B0_GR_omk{o}")
    a0 = L(f"A0_GR_omk{o}")["cl_unlensed_scalar"]
    for m in (0, 1, 2, 3):
        x = L(f"MASK{m}_GR_omk{o}")
        s = stats(x, b0, "cl_unlensed_scalar")
        extra = f"  == A0 bitwise: {np.array_equal(x['cl_unlensed_scalar'], a0)}" if m == 0 else ""
        print(f"omk={o:>6s} {labels[m]:>18s} vs B0: TT(l=2) {s['TT_l2']:.3e} max TT {s['TT_max']:.3e} "
              f"TE/sqrt {s['TE_max']:.3e} EE {s['EE_max']:.3e}{extra}")

"""Verifier (ic1-pidot-eftpid-typo): is the A0 vs FIX P(k) difference for 5_Kmimic_1 stable under accuracy_boost=2?"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

FIX = harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_pidfix"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1typo" / "accboost"
KS = (1.0, 2.0, 5.0, 10.0, 20.0)


def job(args):
    tag, tree, model, over = args
    out = OUT / f"{tag}.npz"
    if not out.exists():
        st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=over, threads=1, timeout=900)
        assert st["status"] == "ok", st
    return tag, harness.load(out)


def pk_at(d, k):
    return float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))


OUT.mkdir(parents=True, exist_ok=True)
jobs = [("GR_ab2", "A0", "1_EFT_GR", {"accuracy_boost": "2"}),
        ("A0_ab2", "A0", "5_Kmimic_1", {"accuracy_boost": "2"}),
        ("FIX_ab2", str(FIX), "5_Kmimic_1", {"accuracy_boost": "2"})]
with ThreadPoolExecutor(3) as ex:
    res = dict(ex.map(job, jobs))
gr = res["GR_ab2"]
for t in ("A0_ab2", "FIX_ab2"):
    print(f"{t}: P/P_GR-1 at k={KS}: " + "  ".join(f"{pk_at(res[t], k) / pk_at(gr, k) - 1:+.4e}" for k in KS))

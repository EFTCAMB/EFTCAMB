"""Verifier (ic1-pidot-eftpid-typo): how much does restoring D2 change each test model?

Runs the given models (default: those missing from t5's pidfix batch plus three overlap models) on FIX
(work/t5-robust/A0_pidfix) with the harness at default settings, then compares with the A0 reference runs
(runs/A0/<model>.npz) and, for the remaining models, uses t5's runs/pidfix outputs. Prints worst Cl / P(k) change.

Usage: python ic1typo_allmodels.py run [models...]   |   python ic1typo_allmodels.py compare
"""

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

FIX = harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_pidfix"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1typo" / "allmodels"
T5 = harness.AUDIT_ROOT / "work" / "t5-robust" / "runs" / "pidfix"
REF = harness.AUDIT_ROOT / "runs" / "A0"
OVERLAP = ["3_RPH_MB_const_1", "2_PEFT_gamma3_const_1", "4_Des_F_R_LCDM_1"]


def default_models():
    have = {p.stem for p in T5.glob("*.npz")}
    return [p.stem for p in harness.model_inis() if p.stem not in have] + OVERLAP


def run(models):
    inis = [harness.PARAMS_DIR / f"{m}.ini" for m in models if not (OUT / f"{m}.npz").exists()]
    st = harness.run_batch(str(FIX), inis, OUT, jobs=4, threads=1, timeout=900)
    for s in st:
        if s["status"] != "ok":
            print("FAILED", s["model_ini"], s["status"])


def pk_rel(a, b):
    k = a["pk_k"]
    kb = b["pk_k"]
    m = (k >= kb[0]) & (k <= kb[-1])
    pb = np.exp(np.interp(np.log(k[m]), np.log(kb), np.log(b["pk"][0])))
    r = np.abs(a["pk"][0][m] / pb - 1)
    i = int(np.argmax(r))
    return float(r[i]), float(k[m][i])


def compare():
    rows = []
    for p in harness.model_inis():
        m = p.stem
        src = OUT / f"{m}.npz"
        tag = "verify"
        if not src.exists():
            src, tag = T5 / f"{m}.npz", "t5"
        if not src.exists() or not (REF / f"{m}.npz").exists():
            rows.append((m, None, None, None, "missing"))
            continue
        a, b = harness.load(src), harness.load(REF / f"{m}.npz")
        c = harness.compare(a, b)
        pkr, pkk = pk_rel(a, b) if "pk" in a and "pk" in b else (np.nan, np.nan)
        rows.append((m, c["worst_cl"], pkr, pkk, tag, c.get("bit_identical")))
    rows.sort(key=lambda r: -(r[2] if r[2] is not None else -1))
    for r in rows:
        if r[1] is None:
            print(f"{r[0]:32s} missing")
            continue
        wc = r[1] if not isinstance(r[1], dict) else max(v for v in r[1].values() if isinstance(v, float))
        print(f"{r[0]:32s} worst_cl {wc if not isinstance(wc, dict) else wc}  worst_pk {r[2]:.3e} at k={r[3]:.3g}"
              f"  src={r[4]} bit_identical={r[5]}")
    (OUT / "compare.json").write_text(json.dumps([[str(x) for x in r] for r in rows], indent=1))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    if sys.argv[1] == "run":
        run(sys.argv[2:] or default_models())
    else:
        compare()

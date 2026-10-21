"""Verify t3-bg-table-resolution-gr-floor: does AccuracyBoost remove the K-mouflage eps2->0 floor? (it should not)."""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness  # noqa: E402

P = harness.PARAMS_DIR
OUT = harness.AUDIT_ROOT / "work" / "verify" / "bgtable_acc"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "2000", "num_redshiftwindows": "0",
        "accuracy_boost": "2", "l_accuracy_boost": "2"}
JOBS = {"gr": (P / "1_EFT_GR.ini", {"omnuh2": "0.0"}),
        "kmof": (P / "5_Kmouflage_1.ini", {"eps2_0": "-1e-08"})}


def job(it):
    n, (ini, ov) = it
    return n, harness.run_model("A0", ini, OUT / f"{n}.npz", overrides={**FAST, **ov}, threads=1, timeout=900)


with ThreadPoolExecutor(2) as ex:
    st = dict(ex.map(job, JOBS.items()))
print({k: v["status"] for k, v in st.items()})
r = harness.compare(OUT / "kmof.npz", OUT / "gr.npz", lmax=1900)
u = r["cl"]["unlensed_scalar"]
print(f"AccuracyBoost=2: dH/H {r['bg_H']:.1e} TT {u['TT']:.1e} EE {u['EE']:.1e} Pk {r['worst_pk']:.1e}")

"""Verifier (spline-integral-convention): is the defect pre-existing? Run the reference and the spline-w designer f(R)
cases on B1 (pre-rebase EFTCAMB on CAMB 1.6.5).
usage: /Users/marco/Software/miniforge3/bin/python verify_spline_integral_b1.py
"""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
from verify_spline_integral_models import CASES, OUT  # noqa: E402


def run(name):
    ini, ov = CASES[name]
    return name, harness.run_model("B1", ini, OUT / f"B1_{name}.npz", overrides=ov, threads=1, timeout=900)


with ThreadPoolExecutor(3) as ex:
    for name, st in ex.map(run, ["DFR_ref", "DFR_s9", "DFR_s10"]):
        errs = [ln.strip()[:120] for ln in st.get("log_tail", []) if "asymptote" in ln or "Error" in ln]
        print(f"B1_{name}: {st['status']} {' / '.join(errs[-2:])}")

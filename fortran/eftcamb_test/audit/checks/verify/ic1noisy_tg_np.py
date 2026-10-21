"""Verifier (ic1-pidot-transient-noisy-coefficients): does the TG_2 single-mode glitch move when the Horndeski
background table resolution changes (model_background_num_points 60000 -> 60001 / 45000)?  A noise-driven IC artefact
should move/change; a deterministic formulation difference should not."""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy" / "tgnp"
JOBS = [(f"tg{m}_np{n}", f"5_hdsk_TG_{m}", {"model_background_num_points": str(n)}) for n in (60001, 45000)
        for m in (1, 2)]


def job(j):
    name, model, ov = j
    out = OUT / f"{name}.npz"
    if not out.exists():
        st = harness.run_model("A0", harness.PARAMS_DIR / f"{model}.ini", out, overrides=ov, threads=1, timeout=900)
        if st["status"] != "ok":
            return name, st["status"]
    return name, harness.load(out)


with ThreadPoolExecutor(4) as ex:
    R = dict(ex.map(job, JOBS))
for n in (60001, 45000):
    a, b = R[f"tg2_np{n}"], R[f"tg1_np{n}"]
    if isinstance(a, str) or isinstance(b, str):
        print(n, a if isinstance(a, str) else "", b if isinstance(b, str) else "")
        continue
    r = a["pk"][0] / b["pk"][0] - 1
    big = np.where(np.abs(r) > 1e-5)[0]
    rt = a["cl_total"][2:2501, 0] / b["cl_total"][2:2501, 0] - 1
    re = a["cl_total"][2:2501, 1] / b["cl_total"][2:2501, 1] - 1
    print(f"np={n}: P(k) TG2/TG1-1: max {np.max(np.abs(r)):.2e}; k with |r|>1e-5: "
          + ", ".join(f"{a['pk_k'][i]:.4g}({r[i]:+.1e})" for i in big)
          + f" | TT max {np.max(np.abs(rt)):.2e} (ell {np.argmax(np.abs(rt)) + 2}), EE max {np.max(np.abs(re)):.2e}"
          f" (ell {np.argmax(np.abs(re)) + 2})")

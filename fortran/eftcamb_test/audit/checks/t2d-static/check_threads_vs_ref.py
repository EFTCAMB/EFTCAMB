"""Thread-race probe, second pass: rerun models on A0 with OMP_NUM_THREADS = N (default 6) and compare with the
2-thread reference runs in audit/runs/A0.  Every k mode is evolved independently and written to its own slot, so all
unlensed / tensor Cls, P(k), sigma8 and the EFT background functions must be bit-identical; only the lensed Cls may
differ at the 1e-15 level (OpenMP reduction order in lensing.f90, upstream CAMB).

Usage: python check_threads_vs_ref.py OUTDIR [--threads N] [model ...]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

REF = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit/runs/A0")
MODELS = ["5_hdsk_TG_2", "5_quint_axion", "5_Kmouflage_1", "5_ADE_1", "5_Horava_eta_1", "5_scg", "2_PEFT_wDE_wCDM_1",
          "4_Des_F_R_LCDM_1", "3_RPH_MB_const_1", "2_PEFT_gamma6_lin_1"]


def main():
    args = sys.argv[1:]
    threads = 6
    if "--threads" in args:
        i = args.index("--threads")
        threads = int(args[i + 1])
        del args[i:i + 2]
    out = Path(args[0])
    models = args[1:] or MODELS
    res = {}
    for m in models:
        st = harness.run_model("A0", harness.PARAMS_DIR / f"{m}.ini", out / f"{m}_t{threads}.npz", threads=threads)
        if st["status"] != "ok":
            res[m] = {"status": st["status"]}
            print(m, res[m], flush=True)
            continue
        a = harness.load(REF / f"{m}.npz")
        b = harness.load(out / f"{m}_t{threads}.npz")
        diffs = {}
        for k in sorted(a):
            if k in b and isinstance(a[k], np.ndarray) and a[k].dtype.kind == "f" and a[k].shape == b[k].shape:
                d = float(np.nanmax(np.abs(a[k] - b[k]))) if a[k].size else 0.0
                if d > 0 or (np.isnan(a[k]) != np.isnan(b[k])).any():
                    diffs[k] = d / max(float(np.nanmax(np.abs(a[k]))), 1e-300)
        nonlens = {k: v for k, v in diffs.items() if k not in ("cl_lensed_scalar", "cl_total")}
        res[m] = {"status": "ok", "diff_fields": diffs, "non_lensing_diff": bool(nonlens)}
        print(m, json.dumps(res[m]), flush=True)
    (out / f"threads_vs_ref_t{threads}.json").write_text(json.dumps(res, indent=1))
    bad = [m for m, r in res.items() if r.get("non_lensing_diff")]
    print(f"models: {len(res)}; with non-lensing differences vs 2-thread reference: {len(bad)} {bad}")


if __name__ == "__main__":
    main()

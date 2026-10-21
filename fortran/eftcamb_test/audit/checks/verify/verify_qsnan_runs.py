"""Harness runs (with get_scale_evolution) for the qs-mu-nan-vanishing-denominator verification.

Usage: python verify_qsnan_runs.py TREE OUTDIR model1 [model2 ...]
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

PAR = Path(__file__).resolve().parents[3] / "parameters"


def main():
    tree, outdir, models = sys.argv[1], Path(sys.argv[2]), sys.argv[3:]
    inis = [PAR / f"{m}.ini" for m in models]
    res = harness.run_batch(tree, inis, outdir, jobs=4, evol_k=[1e-3, 0.01, 0.1, 1.0], threads=1, timeout=900,
                            overrides={"get_tensor_cls": "F"})
    for m, r in zip(models, res):
        print(f"{m:28s} {r['status']:10s} {r['elapsed']}s", r.get("log_tail", [""])[-3:] if r["status"] != "ok" else "")


if __name__ == "__main__":
    main()

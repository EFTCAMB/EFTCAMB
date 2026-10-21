"""EFT-specific numerical knob: background interpolation tables (model_background_num_points) x FACTOR.

python run_nb.py LABEL FACTOR [--jobs 3] [--tree A0] [--extra key=value,...]
Default accuracy, num_redshiftwindows=0 (same as the d1 runs). Only models whose EFT model reads
model_background_num_points are run; the ini value (if any) or the model default is multiplied by FACTOR.
"""

from __future__ import annotations

import argparse
import fnmatch
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t10b-convergence"
# model-file defaults of model_background_num_points by test-ini family
DEFAULTS = {
    "2_PEFT_*": 6000,
    "3_RPH_*": 6000,
    "4_Des_F_R_*": 6000,
    "5_ADE_*": 1500,
    "5_Kmouflage_*": 1000,
    "5_Kmimic_*": 1000,
    "5_quint_*": 1100,
    "5_scg": 10000,
    "5_hdsk_*": 6000,
}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("label")
    p.add_argument("factor", type=float)
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--tree", default="A0")
    p.add_argument("--models", default="*")
    p.add_argument("--extra", default="num_redshiftwindows=0")
    a = p.parse_args()
    extra = dict(s.split("=", 1) for s in a.extra.split(",") if s)
    outdir = WORK / "runs" / a.label
    outdir.mkdir(parents=True, exist_ok=True)
    tasks = []
    for ini in harness.model_inis():
        if not any(fnmatch.fnmatch(ini.stem, pt) for pt in a.models.split(",")):
            continue
        base = next((n for pat, n in DEFAULTS.items() if fnmatch.fnmatch(ini.stem, pat)), None)
        if base is None:
            continue
        base = int(harness.parse_ini(ini).get("model_background_num_points", base))
        if (outdir / f"{ini.stem}.json").exists():
            continue
        ov = dict(extra, model_background_num_points=str(int(base * a.factor)))
        tasks.append((ini, ov))
    print(f"{len(tasks)} tasks", flush=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(harness.run_model, a.tree, ini, outdir / f"{ini.stem}.npz", overrides=ov, threads=2,
                          timeout=5400): (ini, ov) for ini, ov in tasks}
        for f in as_completed(futs):
            ini, ov = futs[f]
            r = f.result()
            print(f"{ini.stem:32s} n={ov['model_background_num_points']:>7s} {r['status']:12s} {r['elapsed']:8.1f}s",
                  flush=True)


if __name__ == "__main__":
    main()

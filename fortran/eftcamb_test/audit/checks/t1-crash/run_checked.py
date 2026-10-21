"""Track t1-crash: run EFTCAMB test models on a runtime-checked build and summarise failures.

Usage:
    python run_checked.py TREE OUTDIR [--pattern '5_*.ini'] [--models a,b] [--jobs 6] [--threads 1] [--timeout 3000]
    python run_checked.py --summary OUTDIR

TREE must be a directory whose basename is AD or ADD so that the harness skips the character getters (F-001).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402


def scan_npz(path: Path) -> dict:
    """Count non-finite values per array in a harness output."""
    bad = {}
    with np.load(path, allow_pickle=False) as f:
        for k in f.files:
            if k == "meta":
                continue
            a = f[k]
            if a.dtype.kind == "f":
                n = int(np.count_nonzero(~np.isfinite(a)))
                if n:
                    bad[k] = n
        meta = json.loads(str(f["meta"])) if "meta" in f.files else {}
    return {"nonfinite": bad, "errors": meta.get("errors", {})}


def backtrace(log: Path) -> list[str]:
    lines = log.read_text(errors="replace").splitlines()
    out = [ln for ln in lines if re.search(r"(Fortran runtime|At line|Program received|#\d+ |Error|error|NaN|Infinity)", ln)]
    return out[:40]


def summary(outdir: Path) -> None:
    for js in sorted(outdir.glob("*.json")):
        info = json.loads(js.read_text())
        name = js.stem
        line = f"{name:40s} {info['status']:16s} {info['elapsed']:8.1f}s"
        npz = js.with_suffix(".npz")
        if info["status"] == "ok" and npz.exists():
            s = scan_npz(npz)
            if s["nonfinite"] or s["errors"]:
                line += f"  NONFINITE={s['nonfinite']} ERR={s['errors']}"
        print(line)
        if info["status"] != "ok":
            for ln in backtrace(Path(info["log"]))[:12]:
                print("      ", ln[:200])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("tree", nargs="?")
    p.add_argument("outdir")
    p.add_argument("--pattern", default="*.ini")
    p.add_argument("--models", default="")
    p.add_argument("--inis", default="", help="explicit comma separated ini paths")
    p.add_argument("--jobs", type=int, default=6)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--timeout", type=float, default=3000)
    p.add_argument("--set", nargs="*", default=[])
    p.add_argument("--summary", action="store_true")
    a = p.parse_args()
    outdir = Path(a.outdir)
    if a.summary:
        summary(outdir)
        return
    if a.inis:
        inis = [Path(x) for x in a.inis.split(",")]
    else:
        inis = harness.model_inis(a.pattern)
        if a.models:
            want = set(a.models.split(","))
            inis = [i for i in inis if i.stem in want]
    over = dict(s.split("=", 1) for s in a.set)
    res = harness.run_batch(a.tree, inis, outdir, jobs=a.jobs, threads=a.threads, timeout=a.timeout, overrides=over)
    for r in res:
        print(f"{Path(r['model_ini']).stem:40s} {r['status']:16s} {r['elapsed']:8.1f}s")


if __name__ == "__main__":
    main()

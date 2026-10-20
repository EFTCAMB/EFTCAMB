"""Run the audit-patched tree (make_patched_tree.py) over test models with a set of AUDIT_* flags and compare.

Usage:
  python flag_batch.py run  TREE OUTDIR "FLAG1,FLAG2|none" [--pattern GLOB] [--set key=value ...] [--jobs N]
  python flag_batch.py cmp  DIR_A DIR_B [--key cl_unlensed_scalar] [--lmax 2500]

cmp prints, per model, the max fractional difference of TT/EE/TE (TE/sqrt(TT EE)) for the chosen Cl key, tensor Cls
and P(k), plus a bit-identity flag.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402


def run(a):
    flags = [] if a.flags == "none" else [f.strip() for f in a.flags.split(",") if f.strip()]
    env = {f"AUDIT_{f.upper()}": "1" for f in flags}
    overrides = dict(s.split("=", 1) for s in a.set)
    inis = harness.model_inis(a.pattern)
    if a.extra_inis:
        inis += [Path(p) for p in a.extra_inis]
    res = harness.run_batch(a.tree, inis, a.outdir, jobs=a.jobs, overrides=overrides, extra_env=env, threads=2)
    for r in res:
        print(f"{Path(r['model_ini']).stem:40s} {r['status']:16s} {r['elapsed']:7.1f}s")


def cmp(a):
    da, db = Path(a.dir_a), Path(a.dir_b)
    rows = {}
    for fa in sorted(da.glob("*.npz")):
        fb = db / fa.name
        if not fb.exists():
            continue
        A, B = harness.load(fa), harness.load(fb)
        c = harness.compare(A, B, lmax=a.lmax)
        row = {"bit": c["bit_identical"], "worst_cl": c["worst_cl"], "worst_pk": c["worst_pk"]}
        row["scal"] = c["cl"].get(a.key.replace("cl_", ""), {})
        row["tens"] = c["cl"].get("tensor", {})
        rows[fa.stem] = row
        s = row["scal"]
        t = row["tens"]
        print(f"{fa.stem:34s} bit={int(row['bit'])} TT={s.get('TT', 0):.2e} EE={s.get('EE', 0):.2e} "
              f"TE={s.get('TE', 0):.2e} tensTT={t.get('TT', 0):.2e} tensBB={t.get('BB', 0):.2e} "
              f"pk={row['worst_pk']:.2e}")
    if a.json:
        Path(a.json).write_text(json.dumps(rows, indent=1))


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("tree")
    r.add_argument("outdir")
    r.add_argument("flags")
    r.add_argument("--pattern", default="*.ini")
    r.add_argument("--extra-inis", nargs="*", default=[])
    r.add_argument("--set", nargs="*", default=[])
    r.add_argument("--jobs", type=int, default=4)
    c = sub.add_parser("cmp")
    c.add_argument("dir_a")
    c.add_argument("dir_b")
    c.add_argument("--key", default="cl_unlensed_scalar")
    c.add_argument("--lmax", type=int, default=2500)
    c.add_argument("--json", default="")
    a = p.parse_args()
    run(a) if a.cmd == "run" else cmp(a)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()

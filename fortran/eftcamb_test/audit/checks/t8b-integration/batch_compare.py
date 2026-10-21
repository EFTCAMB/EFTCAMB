"""t8b: run all test models on a patched tree copy and compare with the A0 reference runs.
Usage: python batch_compare.py TREE_PATH OUTDIR [pattern] [jobs]"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

REF = harness.AUDIT_ROOT / "runs" / "A0"

tree, outdir = sys.argv[1], Path(sys.argv[2])
pattern = sys.argv[3] if len(sys.argv) > 3 else "*.ini"
jobs = int(sys.argv[4]) if len(sys.argv) > 4 else 4
inis = harness.model_inis(pattern)
todo = [i for i in inis if not (outdir / f"{i.stem}.npz").exists()]
harness.run_batch(tree, todo, outdir, jobs=jobs)
summary = {}
for ini in inis:
    a, b = outdir / f"{ini.stem}.npz", REF / f"{ini.stem}.npz"
    st = json.loads((outdir / f"{ini.stem}.json").read_text())["status"] if (outdir / f"{ini.stem}.json").exists() else "?"
    if a.exists() and b.exists():
        c = harness.compare(a, b)
        summary[ini.stem] = {"status": st, "worst_cl": c["worst_cl"], "worst_pk": c["worst_pk"],
                             "bg_H": c.get("bg_H"), "bit_identical": c["bit_identical"],
                             "cl_total": c["cl"].get("total"), "sigma8": c.get("sigma8")}
    else:
        summary[ini.stem] = {"status": st, "missing": True}
(outdir / "summary.json").write_text(json.dumps(summary, indent=1))
for k, v in summary.items():
    if v.get("missing"):
        print(f"{k:32s} {v['status']} MISSING")
    else:
        print(f"{k:32s} {v['status']:6s} cl={v['worst_cl']:.2e} pk={v['worst_pk']:.2e} H={v['bg_H']:.1e} bit={v['bit_identical']}")

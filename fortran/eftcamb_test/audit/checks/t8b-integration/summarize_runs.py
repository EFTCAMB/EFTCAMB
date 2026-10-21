"""t8b: summarise finished runs of a patched-tree batch against the A0 reference runs (no new CAMB runs).
Usage: python summarize_runs.py OUTDIR"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

REF = harness.AUDIT_ROOT / "runs" / "A0"
outdir = Path(sys.argv[1])
rows = []
for npz in sorted(outdir.glob("*.npz")):
    ref = REF / npz.name
    st = json.loads(npz.with_suffix(".json").read_text())["status"] if npz.with_suffix(".json").exists() else "?"
    if not ref.exists():
        continue
    c = harness.compare(npz, ref)
    tot = c["cl"].get("total", {})
    rows.append((npz.stem, st, c["worst_cl"], c["worst_pk"], c.get("sigma8"), tot))
rows.sort(key=lambda r: -max(r[2], r[3]))
for name, st, wcl, wpk, s8, tot in rows:
    s8s = f"{s8:.1e}" if isinstance(s8, float) else str(s8)
    print(f"{name:34s} {st:4s} worst_cl={wcl:.2e} worst_pk={wpk:.2e} sigma8={s8s} "
          + " ".join(f"{k}={v:.1e}" for k, v in tot.items()))
print("n =", len(rows))

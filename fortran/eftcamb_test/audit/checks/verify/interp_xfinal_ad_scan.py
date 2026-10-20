"""Scan every shipped model ini on the bounds-checked tree AD and classify where each run stops.

AD still contains F-002 (results.f90:2168), so a run that passes EFTCAMB init aborts in thermo init; a run that hits the
equispaced-interpolation overflow aborts earlier at 02_equispaced_interpolation_linear_1D.f90.
Usage: python interp_xfinal_ad_scan.py [start end]   (resumable: skips models whose .json exists)
"""

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness  # noqa: E402

OUT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/interp_xfinal/ad_scan")
OUT.mkdir(parents=True, exist_ok=True)
inis = harness.model_inis()
a, b = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (0, len(inis))
todo = [p for p in inis[a:b] if not (OUT / f"{p.stem}.json").exists()]
if todo:
    harness.run_batch("AD", todo, OUT, jobs=4, threads=1, timeout=480)


def classify(log: str) -> str:
    m = re.search(r"At line (\d+) of file \S*/([\w.]+\.f90)", log)
    if m:
        return f"{m.group(2)}:{m.group(1)}"
    if "Fortran runtime error" in log or "Floating-point exception" in log or "SIGFPE" in log:
        return "other-runtime-error"
    return "no-abort"


rows = []
for p in inis:
    js = OUT / f"{p.stem}.json"
    if not js.exists():
        continue
    info = json.loads(js.read_text())
    log = Path(info["log"]).read_text(errors="replace")
    where = classify(log)
    bt = ""
    if "02_equispaced" in where:
        frames = re.findall(r"#\d+\s+0x[0-9a-f]+ in (\S+)", log)
        bt = " <- ".join(f for f in frames[:8])
    rows.append((p.stem, info["status"], where, bt))
for r in rows:
    print(f"{r[0]:40s} {r[1]:14s} {r[2]:50s} {r[3]}")
hits = [r[0] for r in rows if "02_equispaced" in r[2]]
print(f"\nscanned {len(rows)}/{len(inis)}; interp-overflow aborts: {len(hits)}: {hits}")

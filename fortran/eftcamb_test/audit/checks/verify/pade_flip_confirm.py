"""Confirm stability-verdict flips (A0 vs one-line Pade fix) in fresh full harness runs, one process per run."""
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
from pade_stab_scan import CONFIGS, INI, OUT, TREES, point_overrides  # noqa: E402

PTS = [(0.0, 0.0, 0.2, -0.2, -0.6), (0.0, 0.0, 0.2, 0.2, 2.0)]
if len(sys.argv) > 1:
    PTS = [tuple(json.loads(s)) for s in sys.argv[1:]]
jobs = []
for i, p in enumerate(PTS):
    for t in TREES:
        ov = dict(CONFIGS["mink"])
        ov.update(point_overrides(p))
        ov["EFTCAMB_feedback_level"] = "2"
        jobs.append((i, p, t, ov))


def one(j):
    i, p, t, ov = j
    s = harness.run_model(TREES[t], INI, OUT / "confirm" / f"pt{i}_{t}.npz", overrides=ov, threads=1)
    log = Path(s["log"]).read_text(errors="replace")
    reasons = sorted({ln.strip()[:90] for ln in log.splitlines() if "instab" in ln.lower() or "Minkowski" in ln
                      or "unstable" in ln.lower()})
    return {"point": p, "tree": t, "status": s["status"], "reasons": reasons[:4]}


with ThreadPoolExecutor(4) as ex:
    for r in ex.map(one, jobs):
        print(json.dumps(r))

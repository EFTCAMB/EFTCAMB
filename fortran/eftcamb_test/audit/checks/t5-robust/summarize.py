"""Track 5 summary: per variation, compare against defaults with the accuracy_boost=2 noise floor.

For each model with an AB2 run, noise = metrics(AB2, default).  A variation result is flagged when its change in any of
(worst scalar Cl, P(k), sigma8) exceeds max(1e-4, noise of that observable).
Prints a table and writes summary.json in the work directory.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import analyse  # noqa: E402
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t5-robust"
VARS = ["pto_x10", "pto_d10", "grt_x10", "grt_d10", "bto_x10", "bto_d10", "stt_x10", "stt_d10", "ic2", "ic3", "pidfix"]
OBS = ("worst_scalar_cl", "pk", "sigma8")


def main() -> None:
    noise = analyse.analyse("AB2") if (WORK / "runs" / "AB2").exists() else {}
    out = {"noise": noise, "vars": {}}
    print("noise floor (AB2 vs AB1):")
    for m, e in noise.items():
        if e["status"] == "ok":
            print(f"  {m:32s} " + " ".join(f"{o}={e.get(o, float('nan')):.2e}" for o in OBS))
        else:
            print(f"  {m:32s} {e['status']}")
    for var in VARS:
        if not (WORK / "runs" / var).exists():
            continue
        r = analyse.analyse(var)
        out["vars"][var] = r
        nflag = 0
        lines = []
        for m, e in r.items():
            if e["status"] != "ok":
                lines.append(f"  {m:32s} {e['status']}")
                nflag += 1
                continue
            nz = noise.get(m, {})
            flags = []
            for o in OBS:
                v = e.get(o)
                thr = max(1e-4, nz.get(o, 0.0) if nz.get("status") == "ok" else 0.0)
                if v is None:
                    continue
                if not (v <= thr):
                    flags.append(f"{o}={v:.2e}>{thr:.1e}")
            e["flags"] = flags
            if flags:
                nflag += 1
            lines.append(f"  {m:32s} {'FLAG ' + ' '.join(flags) if flags else 'ok'}")
        print(f"=== {var}: {len(r)} models, {nflag} flagged/failed")
        print("\n".join(lines))
    (WORK / "summary.json").write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()

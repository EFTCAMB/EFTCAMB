"""Track 2c: summarise the 300-iteration memory runs written by t2c_memory.py (no CAMB runs).

samples = [iteration, rss_bytes, time_s, malloc_in_use_bytes, malloc_blocks]; growth per iteration is a least-squares
slope over the second half of the run (after allocator / Python warm-up), plus the total change first->last sample.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

M = harness.AUDIT_ROOT / "work" / "t2-state" / "memory"


def slope(it: np.ndarray, y: np.ndarray) -> float:
    h = it >= it.max() / 2
    return float(np.polyfit(it[h], y[h], 1)[0]) if h.sum() >= 3 else float("nan")


def main() -> None:
    print(f"{'run':34s} {'n':>4s} {'rss0':>7s} {'rssN':>7s} {'dRSS/it 2nd half':>17s} {'dMalloc/it 2nd half':>20s} "
          f"{'malloc total':>13s} {'blocks/it':>10s}")
    for p in sorted(M.glob("*/*_300.json")):
        d = json.loads(p.read_text())
        s = np.array(d["samples"], dtype=float)
        if len(s) < 5:
            continue
        it, rss = s[:, 0], s[:, 1] / 1024
        line = f"{p.parent.name + '/' + p.stem:34s} {int(it[-1]) + 1:4d} {rss[0] / 1024:6.0f}M {rss[-1] / 1024:6.0f}M " \
               f"{slope(it, rss):14.1f}KiB"
        if s.shape[1] >= 5:
            mal, blk = s[:, 3] / 1024, s[:, 4]
            line += f" {slope(it, mal):17.2f}KiB {(mal[-1] - mal[0]) / 1024:11.1f}M {slope(it, blk):10.0f}"
        print(line + (f"  errors={d.get('errors')}" if d.get("errors") else ""))


if __name__ == "__main__":
    main()

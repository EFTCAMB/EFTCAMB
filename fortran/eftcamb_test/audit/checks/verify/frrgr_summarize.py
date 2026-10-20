"""Summarise existing frrgr threshold-scan outputs: P_fR/P_GR - 1 at k = 0.31, 1, 3.1, 7.3 and RGR time from logs.

Usage: python frrgr_summarize.py [dir]
"""

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

D = Path(sys.argv[1]) if len(sys.argv) > 1 else harness.AUDIT_ROOT / "work" / "verify" / "frrgr"


def main() -> None:
    for tree in ("A0", "B1"):
        gr = D / f"{tree}_GR.npz"
        if not gr.exists():
            continue
        ref = harness.load(gr)
        k = ref["pk_k"]
        ik = [int(np.argmin(abs(k - x))) for x in (0.31, 1.0, 3.1, 7.3)]
        print(f"[{tree}]  k = " + " ".join(f"{k[i]:.3g}" for i in ik))
        for f in sorted(D.glob(f"{tree}_B0*.npz")):
            d = harness.load(f)
            rgr = re.findall(r"Return to GR time: *([-+0-9.E]+)", f.with_suffix(".log").read_text())
            r = d["pk"][0] / ref["pk"][0] - 1
            print(f"  {f.stem:28s} RGR={float(rgr[-1]) if rgr else float('nan'):7.4f} "
                  + " ".join(f"{r[i]:+11.3e}" for i in ik))


if __name__ == "__main__":
    main()

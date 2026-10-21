"""t3-gr (b): summarise b1_report json: GR-point deviations and approach-to-GR smoothness.

python b2_analyze.py [report.json]
For each family prints worst deviations (Cl / P(k) / sigma8 / H) of the GR point, the forced-EFT-path GR point and
each scan value versus B0 (plain CAMB) and A0g (code under audit at EFTflag=0, IntTolBoost=2), and flags
  - GR point deviation > 1e-4
  - non-smooth approach: deviation at a coupling below 1e-6 larger than 10x (GR-point deviation + 1e-6)
"""

import json
import sys
from pathlib import Path

from common import WORK

THR = 1e-4


def worst(row, ref):
    d = row.get(ref)
    if not d:
        return None
    vals = [d.get(k) for k in ("worst_cl", "pk", "s8", "H")]
    return max(v for v in vals if v is not None)


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else WORK / "b1_report_all.json"
    rep = json.loads(path.read_text())
    flags = []
    for fam, rows in rep.items():
        print(f"== {fam}")
        for tag, row in rows.items():
            if row.get("status") != "ok":
                print(f"   {tag:10s} {row.get('status')} {row.get('log_tail', '')}")
                continue
            parts = []
            for ref in ("B0ref", "A0g", "B0w"):
                if ref in row:
                    d = row[ref]
                    parts.append(f"{ref}: cl={d['worst_cl']:.1e} (TT={d['TT'] or 0:.1e} EE={d['EE'] or 0:.1e} "
                                 f"BB={d['BB'] or 0:.1e} PP={d['PP'] or 0:.1e}) pk={d['pk']:.1e} "
                                 f"s8={d['s8'] or 0:.1e} H={d['H'] or 0:.1e}")
            print(f"   {tag:10s} RGR={row.get('RGR_time')}  " + " | ".join(parts))
        gr = rows.get("GR", {})
        if gr.get("status") == "ok":
            for ref in ("B0ref", "A0g"):
                w = worst(gr, ref)
                if w is not None and w > THR:
                    flags.append((fam, "GR", ref, w))
        g0 = worst(gr, "A0g") if gr.get("status") == "ok" else None
        for tag, row in rows.items():
            if not tag.startswith("s") or row.get("status") != "ok":
                continue
            v = float(tag[1:].replace("m", "-").replace("p", "+"))
            w = worst(row, "A0g")
            if g0 is not None and w is not None and 0 < v <= 1e-6 and w > 10 * (g0 + 1e-6):
                flags.append((fam, tag, "A0g-jump", w))
    print("\nFLAGS (deviation > 1e-4 at GR point, or non-smooth approach):")
    for f in flags:
        print("  ", f)


if __name__ == "__main__":
    main()

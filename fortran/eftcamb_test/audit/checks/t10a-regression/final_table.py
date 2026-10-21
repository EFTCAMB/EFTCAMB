"""Per-model table of the worst inherited-corrected A0-B1 deviation (double ratio vs the matching GR reference).

Reads the JSON written by double_ratio.py and, optionally, the convergence JSON of convergence.py, and prints one
line per model: TT, EE, TE, BB(lensed), tBB(tensor), PP, P(k) and a classification:
  ok        all <= 1e-4
  gr-limit  GR-equivalent model, all <= 1e-8
  known-nc  known not converged at default accuracy (brief): 3_RPH_MB_const_1, 5_Kmimic_2, 5_Horava_SolSyst_*
  >1e-4     needs the convergence classification (see convergence.py output)

Usage: python final_table.py DR_JSON [--conv CONV_JSON] [--json OUT]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

COLS = ["TT", "EE", "TE", "BB_lens", "tBB", "PP", "Pk"]
KNOWN_NC = ("3_RPH_MB_const_1", "5_Kmimic_2", "5_Horava_SolSyst_1", "5_Horava_SolSyst_2", "5_Horava_SolSyst_3")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("dr")
    p.add_argument("--conv", default="")
    p.add_argument("--json", default="")
    a = p.parse_args()
    d = json.loads(Path(a.dr).read_text())["models"]
    conv = json.loads(Path(a.conv).read_text()) if a.conv else {}
    rows = {}
    print(f"{'model':30s} " + " ".join(f"{c:>8s}" for c in COLS) + "  gr_ref          class")
    for m, rec in d.items():
        if "double" not in rec or m == "1_EFT_GR":
            continue
        dd = rec["double"]
        vals = {c: (dd[c]["max"] if isinstance(dd.get(c), dict) and "max" in dd[c] else float("nan")) for c in COLS}
        worst = max(v for v in vals.values() if v == v)
        if m in KNOWN_NC:
            cls = "known-nc"
        elif worst <= 1e-8:
            cls = "gr-limit"
        elif worst <= 1e-4:
            cls = "ok"
        else:
            cls = ">1e-4"
            if m in conv:
                c = conv[m]
                hi = max(v["D_hi"] for v in c.values() if isinstance(v, dict))
                cls += f" D_ab2={hi:.1e}"
        rows[m] = {"vals": vals, "gr_ref": rec.get("gr_ref"), "class": cls}
        print(f"{m:30s} " + " ".join(f"{vals[c]:8.1e}" for c in COLS) + f"  {rec.get('gr_ref', ''):15s} {cls}")
    if a.json:
        Path(a.json).write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()

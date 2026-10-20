"""Track 5: size of the EFT deviations from GR around the perturbation turn-on time of every test model.

For each model reads the default A0 run (eft_* arrays on a = logspace(-8, 0, 401)) and the RGR time printed in the run
log, and reports max over the dimensionless EFT functions (Omega, alpha_B, alpha_M, alpha_T, Gamma_1..6 with
Gamma_2 in units of H0) and the dark-energy density fraction |grhoq| / grho_matter at a = turn-on / 10, turn-on,
turn-on * 10.

Usage: python turnon_magnitudes.py [--json out.json]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

REF = harness.AUDIT_ROOT / "runs" / "A0"
FIELDS = ["EFTOmegaV", "alphaB", "alphaM", "alphaT", "EFTGamma1V", "EFTGamma2V", "EFTGamma3V", "EFTGamma4V",
          "EFTGamma5V", "EFTGamma6V"]


def rgr_time(model: str) -> float | None:
    log = (REF / f"{model}.log").read_text(errors="replace")
    m = re.findall(r"Return to GR time:\s*([0-9.Ee+-]+)", log)
    return float(m[-1]) if m else None


def at(d: dict, field: str, a: float) -> float:
    x = d["eft_a_grid"]
    y = d.get(f"eft_{field}")
    if y is None:
        return 0.0
    return float(np.interp(np.log(a), np.log(x), np.abs(y)))


def deviation(d: dict, a: float) -> dict:
    out = {f: at(d, f, a) for f in FIELDS}
    rho = sum(at(d, f, a) for f in ("grhob_t", "grhoc_t", "grhor_t", "grhog_t", "grhonu_tot"))
    out["DEfrac"] = at(d, "grhoq", a) / rho if rho > 0 else 0.0
    out["max_fn"] = max(out[f] for f in FIELDS)
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--json", default="")
    a = p.parse_args()
    res = {}
    for ini in harness.model_inis():
        m = ini.stem
        d = harness.load(REF / f"{m}.npz")
        keys = harness.parse_ini(ini)
        pto = float(keys.get("EFTCAMB_turn_on_time", "1e-2"))
        r = rgr_time(m) if d["meta"].get("EFTflag", 0) else None
        entry = {"pto": pto, "rgr": r}
        if "eft_a_grid" in d:
            for lab, aa in (("lo", pto / 10), ("on", r if r and r < 1 else pto), ("hi", min(pto * 10, 1.0))):
                entry[lab] = deviation(d, aa)
        res[m] = entry
        if "on" in entry:
            print(f"{m:32s} pto={pto:.0e} rgr={r if r else 0:.3e}  maxfn(lo,on,hi)="
                  f"{entry['lo']['max_fn']:.1e},{entry['on']['max_fn']:.1e},{entry['hi']['max_fn']:.1e}  "
                  f"DEfrac={entry['lo']['DEfrac']:.1e},{entry['on']['DEfrac']:.1e},{entry['hi']['DEfrac']:.1e}")
        else:
            print(f"{m:32s} (GR)")
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

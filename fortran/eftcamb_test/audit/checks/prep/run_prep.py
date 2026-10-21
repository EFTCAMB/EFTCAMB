"""Prep agent: run candidate inis on the A0 tree and summarise (status, parameters read, deviation from GR).

Usage:
    python run_prep.py INI [INI ...] [--jobs 3] [--outdir DIR]

For every ini: harness.run_model on A0 (threads=2), then prints status, model name, parameters read back from
Fortran, max |EFT function| over the a-grid (how far from GR), and the worst Cl / P(k) deviation from the GR
reference run (runs/A0/1_EFT_GR.npz).
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

GR_REF = harness.AUDIT_ROOT / "runs" / "A0" / "1_EFT_GR.npz"
DEFAULT_OUT = harness.AUDIT_ROOT / "work" / "prep" / "runs"
EFT_FIELDS = ("EFTOmegaV", "EFTGamma1V", "EFTGamma2V", "EFTGamma3V", "EFTGamma4V", "EFTGamma5V", "EFTGamma6V")


def summarise(ini: Path, outdir: Path, evol_k: list[float] | None = None) -> dict:
    out = outdir / f"{ini.stem}.npz"
    st = harness.run_model("A0", ini, out, threads=2, evol_k=evol_k)
    res = {"ini": str(ini), "status": st["status"], "elapsed": st["elapsed"], "npz": str(out)}
    if st["status"] != "ok":
        res["log_tail"] = st.get("log_tail", [])[-15:]
        return res
    d = harness.load(out)
    meta = d["meta"]
    res["model_name"] = meta.get("eft_model_name")
    res["errors"] = meta.get("errors")
    res["eft_params"] = {k: v for k, v in (meta.get("eft_params") or {}).items()
                         if k not in ("output_root",) and not k.startswith("EFTCAMB_")}
    res["eft_max"] = {}
    for f in EFT_FIELDS:
        key = f"eft_{f}"
        if key in d:
            v = d[key]
            res["eft_max"][f] = float(np.nanmax(np.abs(v))) if np.any(np.isfinite(v)) else float("nan")
            res["eft_max"][f + "_nan"] = int(np.sum(~np.isfinite(v)))
    if GR_REF.exists():
        c = harness.compare(d, harness.load(GR_REF))
        res["vs_GR"] = {"worst_cl": c["worst_cl"], "worst_pk": c["worst_pk"],
                        "TT": c["cl"].get("total", {}).get("TT"), "sigma8": c.get("sigma8"), "bg_H": c.get("bg_H")}
    res["sigma8"] = float(d["sigma8"][0]) if "sigma8" in d else None
    return res


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("inis", nargs="+")
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--outdir", default=str(DEFAULT_OUT))
    p.add_argument("--json", default="")
    a = p.parse_args()
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    inis = [Path(x).resolve() for x in a.inis]
    with ThreadPoolExecutor(min(a.jobs, 3)) as ex:
        results = list(ex.map(lambda i: summarise(i, outdir), inis))
    for r in results:
        print(json.dumps(r, indent=1, default=str))
    if a.json:
        Path(a.json).write_text(json.dumps(results, indent=1, default=str))


if __name__ == "__main__":
    main()

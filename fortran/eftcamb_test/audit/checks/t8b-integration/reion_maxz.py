"""t8b check: effect of the EFTCAMB change of TTanhReionization%max_redshift (CAMB default 50 -> 1000).

Runs the GR ini (EFTflag=0) and returns the reionization redshift solved from tau and the CMB spectra,
optionally overriding Reion.max_redshift from Python. Compare A0 (default 1000), A0 with 50, and B0 (upstream).
Usage: python reion_maxz.py   (runs the three cases and prints the comparison)
"""

from __future__ import annotations

import sys

import numpy as np


def main(args):
    import camb

    from t8b_common import read_model

    cp = read_model(camb, args["ini"])
    if args.get("max_redshift") is not None:
        cp.Reion.max_redshift = args["max_redshift"]
    res = camb.get_results(cp)
    cl = res.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"]
    return {
        "max_redshift": cp.Reion.max_redshift,
        "zre": res.Params.Reion.redshift,
        "tau": res.get_derived_params().get("tau", None) if False else cp.Reion.optical_depth,
        "cl": cl[:2501].tolist(),
    }


if __name__ == "__main__":
    import json

    from t8b_common import HERE, WORK, run_in_tree

    ini = sys.argv[1] if len(sys.argv) > 1 else "1_EFT_GR.ini"
    cases = {
        "A0_default": ("A0", None),
        "A0_maxz50": ("A0", 50.0),
        "B0_default": ("B0", None),
        "B0_maxz1000": ("B0", 1000.0),
    }
    out = {k: run_in_tree(t, HERE / "reion_maxz.py", {"ini": ini, "max_redshift": m}) for k, (t, m) in cases.items()}
    cl = {k: np.array(v["cl"]) for k, v in out.items()}

    def rel(a, b):
        r = {}
        for i, lab in enumerate(("TT", "EE", "BB")):
            m = b[2:, i] != 0
            r[lab] = float(np.max(np.abs(a[2:, i][m] / b[2:, i][m] - 1))) if m.any() else 0.0
        norm = np.sqrt(np.abs(b[2:, 0] * b[2:, 1]))
        r["TE"] = float(np.max(np.abs(a[2:, 3] - b[2:, 3]) / norm))
        return r

    summary = {k: {"max_redshift": v["max_redshift"], "zre": v["zre"]} for k, v in out.items()}
    summary["A0_default_vs_B0_default"] = rel(cl["A0_default"], cl["B0_default"])
    summary["A0_maxz50_vs_B0_default"] = rel(cl["A0_maxz50"], cl["B0_default"])
    summary["A0_default_vs_B0_maxz1000"] = rel(cl["A0_default"], cl["B0_maxz1000"])
    summary["A0_maxz50_bit_identical_B0"] = bool(np.array_equal(cl["A0_maxz50"], cl["B0_default"]))
    summary["A0_default_bit_identical_B0_maxz1000"] = bool(np.array_equal(cl["A0_default"], cl["B0_maxz1000"]))
    (WORK / "out" / "reion_maxz.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))

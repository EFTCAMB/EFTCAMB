"""EFTflag=0 (GR) consistency of the code under audit (A0) against plain CAMB v2 (B0), flat and non-flat.

The EFTCAMB polterddot replacement in equations.f90 (derivs, output block) is applied for EFTflag=0 too and drops the
curvature factors Kf(2); in flat space it is algebraically identical to the CAMB expression. This script measures
the A0-B0 difference for omk = 0 and omk != 0 with the reionization max_redshift aligned (50 on both).

Usage: python gr_consistency.py OUTDIR
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
base = harness.PARAMS_DIR / "base_params.ini"
cases = {"flat": {}, "omk+0.05": {"omk": "0.05"}, "omk-0.05": {"omk": "-0.05"}}
res = {}
for name, ov in cases.items():
    ov = dict(ov, EFTflag="0", get_tensor_cls="F", num_redshiftwindows="0")
    ra = harness.run_model("A0", base, out / f"A0_{name}.npz", overrides=ov, python_set={"Reion.max_redshift": 50.0})
    rb = harness.run_model("B0", base, out / f"B0_{name}.npz", overrides=ov)
    if ra["status"] != "ok" or rb["status"] != "ok":
        res[name] = {"A0": ra["status"], "B0": rb["status"]}
        continue
    c = harness.compare(out / f"A0_{name}.npz", out / f"B0_{name}.npz")
    res[name] = {"bit_identical": c["bit_identical"], "worst_cl": c["worst_cl"], "worst_pk": c["worst_pk"],
                 "cl_total": c["cl"].get("total"), "cl_unlensed_scalar": c["cl"].get("unlensed_scalar")}
print(json.dumps(res, indent=1))
(out / "gr_consistency.json").write_text(json.dumps(res, indent=1))

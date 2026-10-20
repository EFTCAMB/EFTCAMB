"""t8b check: results.f90 ThermoData_Init doubles the background sampling (ninverse) only when
CP%DarkEnergy%is_cosmological_constant is false. With EFTflag /= 0 CP%DarkEnergy is the never-initialised
placeholder (w=-1), so EFT models with non-LCDM backgrounds never get the doubling. For the models most sensitive
to it (from the A0_ninv batch), compare A0 (no doubling) and A0_ninv (doubling) against a boosted reference
(A0, Accuracy.BackgroundTimeStepBoost = 4). Usage: python ninv_reference.py"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "t8b-integration"
models = sys.argv[1:] or ["5_Kmimic_1", "5_Kmimic_3", "5_ADE_1"]
outdir = W / "runs_btsb4"
outdir.mkdir(parents=True, exist_ok=True)
for m in models:
    npz = outdir / f"{m}.npz"
    if not npz.exists():
        st = harness.run_model("A0", harness.PARAMS_DIR / f"{m}.ini" if hasattr(harness, "PARAMS_DIR") else
                               Path(__file__).resolve().parents[3] / "parameters" / f"{m}.ini", npz,
                               python_set={"Accuracy.BackgroundTimeStepBoost": 4.0})
        print(m, st["status"])
res = {}
for m in models:
    ref = outdir / f"{m}.npz"
    a0 = harness.AUDIT_ROOT / "runs" / "A0" / f"{m}.npz"
    ni = W / "runs_ninv" / f"{m}.npz"
    r = {}
    for name, f in (("A0_default", a0), ("A0_ninv_doubled", ni)):
        if f.exists() and ref.exists():
            c = harness.compare(f, ref)
            r[name] = {"worst_cl": c["worst_cl"], "worst_pk": c["worst_pk"], "sigma8": c.get("sigma8")}
    res[m] = r
    print(m, json.dumps(r))
(W / "out" / "ninv_reference.json").write_text(json.dumps(res, indent=1))

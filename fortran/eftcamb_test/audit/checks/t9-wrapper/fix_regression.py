"""Track 9: regression check of the candidate-fix tree (work/t9-wrapper/A0fix) against the A0 reference runs.

Runs a few models through harness.run_model on A0fix and compares with runs/A0/<model>.npz (bit identity expected:
the candidate fixes only touch the python interface, error handling and a guard in CAMBdata_SetParams).

Usage: python fix_regression.py [model ...]
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

from t9lib import AUDIT, PARAMS, WORK  # noqa: E402

TREE = WORK / "A0fix"
OUT = WORK / "fixreg"
MODELS = sys.argv[1:] or ["2_PEFT_Omega_power_1", "4_Des_F_R_LCDM_2", "3_RPH_KT_power_1", "1_EFT_GR"]

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    res = harness.run_batch(str(TREE), [PARAMS / f"{m}.ini" for m in MODELS], str(OUT), jobs=2)
    for r in res:
        m = Path(r["model_ini"]).stem
        if r["status"] != "ok":
            print(f"{m:28s} {r['status']}")
            continue
        c = harness.compare(harness.load(OUT / f"{m}.npz"), harness.load(AUDIT / "runs" / "A0" / f"{m}.npz"))
        print(f"{m:28s} bit_identical={c['bit_identical']} worst_cl={c['worst_cl']:.2e} worst_pk={c['worst_pk']:.2e}"
              f" missing={c.get('missing')}")
    print(json.dumps({"tree": str(TREE)}))

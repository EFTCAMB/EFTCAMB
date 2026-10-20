"""Independent analysis of IC2/IC3 runs (any directory of harness npz) against the IC1 A0 reference runs."""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

REF = harness.AUDIT_ROOT / "runs" / "A0"
for d in sys.argv[1:]:
    d = Path(d)
    print(f"== {d}")
    for js in sorted(d.glob("*.json")):
        st = json.loads(js.read_text())
        model = js.stem
        ini = harness.parse_ini(harness.PARAMS_DIR / f"{model}.ini")
        ton = ini.get("EFTCAMB_turn_on_time", "0.01(default)")
        row = f"{model:28s} ton={ton:>8s} status={st['status']:12s}"
        if st["status"] == "ok":
            a = harness.load(js.with_suffix(".npz"))
            b = harness.load(REF / f"{model}.npz")
            cl = a["cl_total"][2:2501]
            nanfrac = float(np.mean(~np.isfinite(cl)))
            if nanfrac > 0:
                row += f" NaN fraction {nanfrac:.2f}"
            else:
                tt = np.max(np.abs(a["cl_total"][2:2501, 0] / b["cl_total"][2:2501, 0] - 1))
                pk = np.max(np.abs(a["pk"][0] / b["pk"][0] - 1))
                row += f" max|dTT/TT|={tt:9.2e} max|dP/P|(z0)={pk:9.2e} sigma8 {float(b['sigma8'][0]):.4f}->{float(a['sigma8'][0]):.4f}"
            errs = a["meta"].get("errors")
            if errs:
                row += f" errors={list(errs)}"
        else:
            tail = " | ".join(x.strip() for x in st.get("log_tail", [])[-6:] if "Error" in x or "error" in x)
            row += " " + tail[:150]
        print(row)

"""Verify ic3-nan-when-omega-gamma3-zero via harness: IC 1/2/3 on minimally-coupled models, IC3 on a nonzero-Omega model."""
import sys, json
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness
OUT = harness.AUDIT_ROOT / "work/verify/ic3nan"
P = harness.REPO / "fortran/eftcamb_test/parameters"
cases = [(m, ic) for m in sys.argv[1].split(",") for ic in map(int, sys.argv[2].split(","))]
for m, ic in cases:
    out = OUT / f"{m}_ic{ic}.npz"
    st = harness.run_model(harness.AUDIT_ROOT / "trees/A0", P / f"{m}.ini", out, overrides={"EFT_IC_type": ic}, threads=1)
    r = {"model": m, "ic": ic, "status": st["status"]}
    if st["status"] == "ok":
        d = harness.load(out)
        cl = d["cl_total"][2:2501]
        r.update(nan_frac=float(np.isnan(cl).mean()), sigma8=float(np.ravel(d["sigma8"])[-1]),
                 TT100=float(cl[98, 0]), errors=(d["meta"].get("errors") if isinstance(d["meta"], dict) else str(d["meta"])[:200]),
                 Omega_a1e8=float(d["eft_EFTOmegaV"][0]) if "eft_EFTOmegaV" in d else None)
    else:
        r["tail"] = st.get("log_tail", "")[-300:]
    print(json.dumps(r), flush=True)

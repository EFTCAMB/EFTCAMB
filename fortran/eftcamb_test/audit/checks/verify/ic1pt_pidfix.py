"""Attribution: typo-only fixed tree (t5-robust/A0_pidfix: EFTpiD1 -> EFTpiD2 in EFTpiDfunction, pi' still from dfridr)."""
import json, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1pt"
PF = str(harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_pidfix")
def job(a):
    name, tree, model, ov = a
    out = OUT / f"{name}.npz"
    if not out.exists():
        st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=ov, threads=1, timeout=600)
        if st["status"] != "ok":
            return name, st["status"]
    return name, harness.load(out)
def pk_at(d, k):
    return float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))
jobs = [("A0_GR", "A0", "1_EFT_GR", {}), ("A0_TG1", "A0", "5_hdsk_TG_1", {}), ("pf_TG2", PF, "5_hdsk_TG_2", {})]
for t in ("1e-9", "1.0000e-08", "1.0200e-08", "1e-7"):
    jobs.append((f"pf_km_{t}", PF, "5_Kmimic_1", {"EFTCAMB_turn_on_time": t}))
with ThreadPoolExecutor(4) as ex:
    res = dict(ex.map(job, jobs))
s = {}
for n, d in res.items():
    if n.startswith("pf_km"):
        s[n] = d if isinstance(d, str) else {f"k={k:g}": pk_at(d, k) / pk_at(res["A0_GR"], k) - 1 for k in (0.1, 1, 5, 10)}
a, b = res["A0_TG1"], res["pf_TG2"]
for j, sp in ((0, "TT"), (1, "EE")):
    q = np.abs(b["cl_total"][2:2501, j] / a["cl_total"][2:2501, j] - 1)
    s[f"pf_TG2_vs_A0_TG1_{sp}"] = (float(q.max()), int(np.argmax(q)) + 2)
print(json.dumps(s, indent=1))

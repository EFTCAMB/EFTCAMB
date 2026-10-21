"""Independent check for finding ic1-pidot-transient-noisy-coefficients.

1. Fine turn-on scan of 5_Kmimic_1 on A0 (1e-8 * {1, 1.02, 1.05, 1.1, 1.2}) and on A0_icsw with pi'_ini = 0:
   a smooth physical transient would give smooth dependence; noisy coefficients give erratic jumps.
2. B1 (pre-rebase) Kmimic turn-on 1e-9/1e-8/1e-7 and TG_2 vs TG_1, for attribution.
3. Smoothness of EFTpiC / EFTpiD1 / EFTpiE from get_eft_functions (different code path from the IC routine).
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1pt"
ICSW = str(harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_icsw")
KS = (0.1, 1.0, 5.0, 10.0)


def job(a):
    name, tree, model, ov, env = a
    out = OUT / f"{name}.npz"
    if not out.exists():
        st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=ov, threads=1, timeout=600,
                               extra_env=env)
        if st["status"] != "ok":
            return name, st["status"]
    return name, harness.load(out)


def pk_at(d, k):
    return float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))


def main():
    part = sys.argv[1]
    jobs = [("A0_GR", "A0", "1_EFT_GR", {}, {}), ("B1_GR", "B1", "1_EFT_GR", {}, {})]
    if part == "scan":
        for f in ("1.0", "1.02", "1.05", "1.1", "1.2"):
            t = f"{float(f) * 1e-8:.4e}"
            jobs.append((f"A0_km_{f}", "A0", "5_Kmimic_1", {"EFTCAMB_turn_on_time": t}, {}))
            jobs.append((f"icsw0_km_{f}", ICSW, "5_Kmimic_1", {"EFTCAMB_turn_on_time": t}, {"EFT_T5_PIDOT0": "1"}))
    elif part == "b1":
        for t in ("1e-9", "1e-8", "1e-7"):
            jobs.append((f"B1_km_{t}", "B1", "5_Kmimic_1", {"EFTCAMB_turn_on_time": t}, {}))
        jobs.append(("B1_TG1", "B1", "5_hdsk_TG_1", {}, {}))
        jobs.append(("B1_TG2", "B1", "5_hdsk_TG_2", {}, {}))
        jobs.append(("A0_TG1", "A0", "5_hdsk_TG_1", {}, {}))
        jobs.append(("A0_TG2", "A0", "5_hdsk_TG_2", {}, {}))
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(job, jobs))
    summ = {}
    for name, d in res.items():
        if "GR" in name and "TG" not in name:
            continue
        if isinstance(d, str):
            summ[name] = d
            continue
        if "_km_" in name:
            gr = res["B1_GR" if name.startswith("B1") else "A0_GR"]
            summ[name] = {f"k={k:g}": pk_at(d, k) / pk_at(gr, k) - 1 for k in KS}
    if part == "b1":
        for tag in ("A0", "B1"):
            a, b = res[f"{tag}_TG1"], res[f"{tag}_TG2"]
            if isinstance(a, str) or isinstance(b, str):
                summ[f"{tag}_TG"] = "fail"
                continue
            r = {}
            for j, s in ((0, "TT"), (1, "EE")):
                x, y = b["cl_total"][2:2501, j], a["cl_total"][2:2501, j]
                q = np.abs(x / y - 1)
                r[s] = (float(q.max()), int(np.argmax(q)) + 2)
            summ[f"{tag}_TG2_vs_TG1"] = r
    print(json.dumps(summ, indent=1))
    (OUT / f"summary_{part}.json").write_text(json.dumps(summ, indent=1))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    main()

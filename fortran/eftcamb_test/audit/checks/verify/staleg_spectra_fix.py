"""Spectra impact of ic23-stale-cache-gammas: run models with EFT_IC_type = 3 (or 2) on work/t5-robust/A0_icdbg with the
release code path (EFT_T5_ICFIX=0) and with the candidate fix (EFT_T5_ICFIX=1, full cache at a_eval), 1 thread, and
also the frozen A0 tree (to check A0_icdbg ICFIX=0 == A0).  usage: python staleg_spectra_fix.py IC model [model ...]"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

ic = sys.argv[1]
models = sys.argv[2:]
OUT = harness.AUDIT_ROOT / "work" / "verify" / "staleg" / "spectra"
ICDBG = harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_icdbg"
jobs = []
for m in models:
    for tag, tree, env in (("A0", "A0", {}), ("rel", ICDBG, {"EFT_T5_ICFIX": "0"}), ("fix", ICDBG, {"EFT_T5_ICFIX": "1"})):
        jobs.append((m, tag, tree, env))


def one(j):
    m, tag, tree, env = j
    out = OUT / f"{m}_ic{ic}_{tag}.npz"
    st = harness.run_model(tree, harness.PARAMS_DIR / f"{m}.ini", out, overrides={"EFT_IC_type": ic}, threads=1,
                           extra_env=env, timeout=1500)
    return (m, tag), st


with ThreadPoolExecutor(4) as ex:
    sts = dict(ex.map(one, jobs))

summary = {}
for m in models:
    row = {t: sts[(m, t)]["status"] for t in ("A0", "rel", "fix")}
    d = {}
    for t in ("A0", "rel", "fix"):
        if row[t] == "ok":
            d[t] = harness.load(OUT / f"{m}_ic{ic}_{t}.npz")
            row[f"finite_{t}"] = bool(np.all(np.isfinite(d[t]["cl_total"][2:2501])) and np.all(np.isfinite(d[t]["pk"])))
            row[f"sigma8_{t}"] = float(np.atleast_1d(d[t]["sigma8"])[0])
    if "A0" in d and "rel" in d:
        c = harness.compare(d["rel"], d["A0"])
        row["rel_vs_A0_bitident"] = c["bit_identical"]
        row["rel_vs_A0_worst_cl"] = c["worst_cl"]
    if "rel" in d and "fix" in d:
        c = harness.compare(d["fix"], d["rel"])
        row.update(fix_vs_rel_TT=c["cl"]["total"]["TT"], fix_vs_rel_EE=c["cl"]["total"]["EE"],
                   fix_vs_rel_TE=c["cl"]["total"]["TE"], fix_vs_rel_pk0=c["pk"][0], fix_vs_rel_worst_cl=c["worst_cl"])
    if not all(row[t] == "ok" for t in ("A0", "rel", "fix")):
        for t in ("A0", "rel", "fix"):
            if row[t] != "ok":
                row[f"tail_{t}"] = [x for x in sts[(m, t)].get("log_tail", []) if "rror" in x][-2:]
    summary[m] = row
    print(m, json.dumps(row, default=str))
(OUT / f"summary_ic{ic}_{'_'.join(models)[:80]}.json").write_text(json.dumps(summary, indent=1, default=str))

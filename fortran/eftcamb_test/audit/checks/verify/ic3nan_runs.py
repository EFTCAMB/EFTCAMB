"""Independent check of finding ic3-nan-when-omega-gamma3-zero via harness.run_model (different code path from the
repro: harness worker, get_results + get_scale_evolution), plus controls and attribution.

Runs (threads=1, <=4 concurrent):
  A0 gamma1_const_2 IC3 with evol_k  -> does the NaN appear before the EFT turn-on (a < 0.01)? => equations.f90 initial()
  A0 gamma1_const_2 IC2 / IC1        -> case 2 and GR are finite (control)
  A0 5_quint_exponential, 4_Des_mc5e_CPL_1, 2_PEFT_wDE_wCDM_1, 2_PEFT_gamma6_const_1 IC3   -> breadth claim
  A0 2_PEFT_Omega_const_2 (Om=0.1), 2_PEFT_gamma3_const_2 (G3=0.1) IC3   -> finite (control: nonzero denominator)
  B1 gamma1_const_2 IC3              -> pre-rebase EFTCAMB, same code: regression or pre-existing?
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

AUD = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
OUT = Path(AUD + "/work/verify/ic3nan/runs")
PAR = Path("/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")

jobs = [
    ("A0", "2_PEFT_gamma1_const_2", 3, [1e-3, 0.05]),
    ("A0", "2_PEFT_gamma1_const_2", 2, None),
    ("A0", "2_PEFT_gamma1_const_2", 1, [1e-3, 0.05]),
    ("A0", "5_quint_exponential", 3, None),
    ("A0", "4_Des_mc5e_CPL_1", 3, None),
    ("A0", "2_PEFT_wDE_wCDM_1", 3, None),
    ("A0", "2_PEFT_gamma6_const_1", 3, None),
    ("A0", "2_PEFT_Omega_const_2", 3, None),
    ("A0", "2_PEFT_gamma3_const_2", 3, None),
    ("B1", "2_PEFT_gamma1_const_2", 3, None),
]
if len(sys.argv) > 1:
    sel = set(int(i) for i in sys.argv[1].split(","))
    jobs = [j for i, j in enumerate(jobs) if i in sel]


def one(job):
    tree, model, ic, evk = job
    out = OUT / f"{tree}_{model}_IC{ic}.npz"
    js = out.with_suffix(".json")
    if js.exists() and out.exists() and len(sys.argv) <= 2:
        st = json.loads(js.read_text())  # reuse a finished run (pass a 2nd arg to force rerun)
    else:
        st = harness.run_model(AUD + f"/trees/{tree}", PAR / f"{model}.ini", out,
                               overrides={"EFT_IC_type": str(ic)}, evol_k=evk, threads=1, timeout=900)
    res = {"tree": tree, "model": model, "IC": ic, "status": st["status"], "elapsed": st["elapsed"]}
    if st["status"] == "ok":
        d = harness.load(out)
        cl = d.get("cl_total")
        if cl is not None:
            c = cl[2:2501]
            res["cl_nan_frac"] = float(np.isnan(c).mean())
            res["TT_l100"] = float(c[98, 0])
        res["sigma8"] = float(d["sigma8"][0]) if "sigma8" in d else None
        res["pk_nan_frac"] = float(np.isnan(d["pk"]).mean()) if "pk" in d else None
        res["errors"] = d["meta"].get("errors")
        res["eft_init_fallback"] = d["meta"].get("eft_init_fallback")
        res["EFT_IC_type_meta"] = (d["meta"].get("eft_params") or {}).get("EFT_IC_type")
        if "evol_clxc" in d:
            a = np.logspace(-6, 0, 300)  # harness a_ev grid (npz evol_a is the cache field 'a', shape [k, a])
            ev = {}
            for f in ("clxc", "clxg", "etak", "pi", "z"):
                key = f"evol_{f}"
                if key in d:
                    arr = d[key]
                    pre = a < 0.01
                    ev[f] = {
                        "nan_frac_pre_turn_on": float(np.isnan(arr[:, pre]).mean()),
                        "nan_frac_post_turn_on": float(np.isnan(arr[:, ~pre]).mean()),
                        "first_val_k0": float(arr[0, 0]),
                    }
            res["evol"] = ev
    else:
        res["log_tail"] = st.get("log_tail", [])[-8:]
    return res


with ThreadPoolExecutor(4) as ex:
    results = list(ex.map(one, jobs))
for r in results:
    print(json.dumps(r))
tag = sys.argv[1].replace(",", "_") if len(sys.argv) > 1 else "all"
(OUT / f"summary_{tag}.json").write_text(json.dumps(results, indent=1))

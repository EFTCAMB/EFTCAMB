"""Candidate fix for the pi'' (pidotdot) read-before-write in derivs: tree FX = A0 + refresh of the Einstein factors
(and cached dz) right after the pi equation (see make notes in the t8a findings). Runs a few models on FX and compares
with A0 (code under audit) and with P + AUDIT_FIX_PIDOTDOT=1 (double derivs call in output()).

Usage: python fx_compare.py
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/t8a-equations")
A0RUNS = harness.AUDIT_ROOT / "runs" / "A0"
PRUNS = W / "runs" / "P_pidot"
OUT = W / "runs" / "FX"
OUT.mkdir(parents=True, exist_ok=True)
MODELS = ["5_Horava_eta_4", "5_Horava_eta_3", "2_PEFT_Omega_const_1", "4_Des_F_R_LCDM_1", "3_RPH_KT_power_1",
          "5_hdsk_TG_2", "2_PEFT_gamma3_const_1", "1_EFT_GR"]


def run(m):
    return harness.run_model(W / "FX", harness.PARAMS_DIR / f"{m}.ini", OUT / f"{m}.npz", threads=2)


with ThreadPoolExecutor(4) as ex:
    st = dict(zip(MODELS, ex.map(run, MODELS)))
res = {}
for m in MODELS:
    if st[m]["status"] != "ok":
        res[m] = st[m]["status"]
        continue
    r = {}
    for tag, ref in (("vs_A0", A0RUNS / f"{m}.npz"), ("vs_P_fixpidotdot", PRUNS / f"{m}.npz")):
        if ref.exists():
            c = harness.compare(OUT / f"{m}.npz", ref)
            r[tag] = {"bit": c["bit_identical"], "unlensed_scalar": c["cl"].get("unlensed_scalar"),
                      "worst_pk": c["worst_pk"]}
    res[m] = r
    print(m, json.dumps(r))
(W / "fx_compare.json").write_text(json.dumps(res, indent=1))

"""Verifier: independent ini-path check of EFTCAMB_turn_on_time <= 0 / NaN / > 1 giving silent GR perturbations.

Runs two test inis through harness.run_model (camb.read_ini path, different from the repro's python
initialize_parameters path) with overrides, and compares TT/EE against the default-turn-on run and an EFTflag=0 run.
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
import numpy as np  # noqa: E402

TREE = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0"
OUT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/turn_on_invalid")
OUT.mkdir(parents=True, exist_ok=True)
PAR = Path(__file__).resolve().parents[3] / "parameters"
MODELS = ["2_PEFT_Omega_const_1", "4_Des_F_R_LCDM_1"]
CASES = {"default": {}, "gr": {"EFTflag": 0}, "on0": {"EFTCAMB_turn_on_time": 0.0},
         "onneg": {"EFTCAMB_turn_on_time": -0.1}, "on2": {"EFTCAMB_turn_on_time": 2.0},
         "on1e30": {"EFTCAMB_turn_on_time": 1e-30}}

jobs = []
for m in MODELS:
    for c, ov in CASES.items():
        ov = dict(ov, feedback_level=2, l_max_scalar=1500)
        jobs.append((m, c, ov))


def go(j):
    m, c, ov = j
    st = harness.run_model(TREE, str(PAR / f"{m}.ini"), str(OUT / f"{m}_{c}.npz"), overrides=ov, threads=1, timeout=900)
    return m, c, st


with ThreadPoolExecutor(4) as ex:
    res = list(ex.map(go, jobs))
status = {(m, c): st for m, c, st in res}
for m in MODELS:
    d = harness.load(str(OUT / f"{m}_default.npz"))
    g = harness.load(str(OUT / f"{m}_gr.npz"))

    def mx(a, b, col):
        x, y = a["cl_total"][2:1500, col], b["cl_total"][2:1500, col]
        return np.max(np.abs(x / y - 1))

    print(f"{m}: GR vs default TT {mx(g, d, 0):.3e}  EE {mx(g, d, 1):.3e}")
    for c in CASES:
        if c in ("default", "gr"):
            continue
        st = status[(m, c)]
        if st["status"] != "ok":
            print(f"   {c:8s} status={st['status']} tail={st.get('log_tail', '')[-150:]!r}")
            continue
        r = harness.load(str(OUT / f"{m}_{c}.npz"))
        meta = r.get("meta")
        try:
            meta = json.loads(str(meta))
        except Exception:
            meta = {}
        print(f"   {c:8s} ok  vs default TT {mx(r, d, 0):.3e}  vs GR TT {mx(r, g, 0):.3e} EE {mx(r, g, 1):.3e}  "
              f"errors={meta.get('errors')}")
        log = Path(st["log"]).read_text(errors="replace")
        for line in log.splitlines():
            if "RGR" in line or "return to GR" in line.lower() or "turn on" in line.lower():
                print("        log:", line.strip()[:120])

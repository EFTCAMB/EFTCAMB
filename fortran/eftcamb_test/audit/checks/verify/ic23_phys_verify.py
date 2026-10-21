"""Verifier: EFT_IC_type 2/3 pi series at turn-on. Runs IC1/IC2/IC3 for several turn-on times via harness."""
import sys, json
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
import numpy as np
from concurrent.futures import ThreadPoolExecutor

ROOT = harness.AUDIT_ROOT
TREE = str(ROOT / "trees" / "A0")
OUT = ROOT / "work" / "verify" / "ic23phys"
PAR = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/"
model = sys.argv[1]
ics = [int(x) for x in sys.argv[2].split(",")]
tons = [float(x) for x in sys.argv[3].split(",")]
jobs = []
for ic in ics:
    for t in tons:
        jobs.append((ic, t))

def go(j):
    ic, t = j
    out = OUT / f"{model}_ic{ic}_on{t:g}.npz"
    if out.exists():
        return j, {"status": "ok(cached)"}
    st = harness.run_model(TREE, PAR + model + ".ini", str(out), overrides={"EFT_IC_type": ic, "EFTCAMB_turn_on_time": t},
                           evol_k=[0.01, 0.1, 1.0], threads=1, timeout=1500)
    return j, st

with ThreadPoolExecutor(4) as ex:
    for j, st in ex.map(go, jobs):
        print(j, st.get("status"), (st.get("log_tail") or "")[-200:] if "ok" not in st.get("status") else "")

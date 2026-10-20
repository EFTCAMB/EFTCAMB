"""Verifier for t3-eft-turnon-before-tauminn-lmax: B1 at l_max 2000, and A0 failing via turn-on time at l_max 2500."""
import re, sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
INI = harness.PARAMS_DIR / "5_ADE_1.ini"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "tauminn_turnon"
CASES = [("B1", {"l_max_scalar": "2000"}), ("B1", {"l_max_scalar": "2500"}),
         ("A0", {"l_max_scalar": "2500", "EFTCAMB_turn_on_time": "5.e-11"}),
         ("A0", {"l_max_scalar": "2500", "EFTCAMB_turn_on_time": "2.e-10"}),
         ("A0", {"l_max_scalar": "2000", "EFTCAMB_turn_on_time": "1.2e-10"})]
def job(c):
    tree, ov = c
    ov = dict(ov, get_tensor_cls="F")
    name = tree + "_" + "_".join(f"{k}{v}" for k, v in ov.items())
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, INI, out, overrides=ov, threads=1, timeout=900)
    msg = re.findall(r"(EFTCAMB starts before[^\n]*)", out.with_suffix(".log").read_text(errors="replace"))
    return name, st["status"], msg[:1]
with ThreadPoolExecutor(3) as ex:
    for r in ex.map(job, CASES):
        print(*r)

"""Does the out-of-bounds read ever fault in release? Run 3_RPH_KT_power_1 on A0 (frozen release) and on the
instrumented release tree with model_background_num_points chosen so the 6 x N-double tables are page-sized / large
allocations (where one element past the end can fall on an unmapped page)."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness  # noqa: E402

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
OUT = AUDIT / "work/verify/interp_xfinal/npts"
INI = Path("/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/3_RPH_KT_power_1.ini")
tree = sys.argv[1] if len(sys.argv) > 1 else "A0"
npts = [int(v) for v in sys.argv[2:]] or [16384, 32768, 65536, 131072]
tag = "instr" if "instr" in tree else Path(tree).name
from concurrent.futures import ThreadPoolExecutor


def one(n):
    return n, harness.run_model(tree, INI, OUT / f"{tag}_N{n}.npz", overrides={"model_background_num_points": str(n)},
                                threads=1, timeout=540)


with ThreadPoolExecutor(4) as ex:
    res = list(ex.map(one, npts))
for n, info in res:
    log = Path(info["log"]).read_text(errors="replace").splitlines()
    nhit = sum(ln.startswith("AUDIT_OOB") for ln in log)
    ex1 = next((ln for ln in log if ln.startswith("AUDIT_OOB")), "")
    print(f"{tag} N={n:7d} status={info['status']:14s} oob_hits={nhit} {ex1[:200]}")
    if info["status"] != "ok":
        print("   tail:", " | ".join(info.get("log_tail", [])[-4:]))

"""Adversarial verifier (t3-extgal-gr-point-background): is the hdd_branch vs allfix C_l difference at S=0 due to nu?

S=0, unlensed scalar only (l <= 2500), omnuh2 = 0.01 and omnuh2 = 0: hdd_branch vs allfix.
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work/verify/extgal_gr"
OUT = W / "nuctl"
INI = REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"
BASE = {"ExtendedGalileon_S": "0", "ExtendedGalileon_q": "100", "get_tensor_cls": "F", "do_lensing": "F",
        "l_max_scalar": "2500", "num_redshiftwindows": "0"}
jobs = [(t, nu) for nu in ("0.01", "0.0") for t in ("hdd_branch", "allfix")]


def one(j):
    t, nu = j
    out = OUT / f"{t}_nu{nu}.npz"
    return j, harness.run_model(W / t, INI, out, overrides={**BASE, "omnuh2": nu}, threads=1, timeout=900), out


with ThreadPoolExecutor(4) as ex:
    res = {j: (st, out) for j, st, out in ex.map(one, jobs)}
for nu in ("0.01", "0.0"):
    a, b = harness.load(res[("hdd_branch", nu)][1]), harness.load(res[("allfix", nu)][1])
    c = harness.compare(a, b)
    print(f"omnuh2={nu}: hdd_branch vs allfix unlensed TT={c['cl']['unlensed_scalar']['TT']:.2e} "
          f"EE={c['cl']['unlensed_scalar']['EE']:.2e} worst_pk={c['worst_pk']:.2e} bg_H={c['bg_H']:.2e}")

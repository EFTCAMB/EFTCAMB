"""Verify t3-stability-dtauda-without-matter (part D): value of grhoa2 filled in derivs (equations.f90:2792) and in
the python bridge (11_EFTCAMB_python.f90:418), as returned by get_scale_evolution / get_eft_functions, vs the
convention used by every consumer (equations.f90 dtauda: grhoa2 = grho_no_de(a) = grhom_t * a^2).

Run: python verify_stab_dtauda_derivs_grhoa2.py
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REPO = HERE.parents[4]
OUT = harness.AUDIT_ROOT / "work" / "verify" / "stabdtauda" / "derivs"
INI = REPO / "fortran" / "eftcamb_test" / "parameters" / "5_Kmouflage_1.ini"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "600", "get_transfer": "F",
        "num_redshiftwindows": "0"}
for tag, tree in (("A0", "A0"), ("R5", str(harness.AUDIT_ROOT / "work" / "t3-gr" / "R5"))):
    out = OUT / f"{tag}_Kmou_nom_evol.npz"
    st = harness.run_model(tree, INI, out, overrides=FAST, evol_k=[0.1], threads=1, timeout=300)
    d = harness.load(out)
    a = d["evol_a"]
    ev_ratio = d["evol_grhoa2"][0] / (d["evol_grhom_t"][0] * a**2)
    ea = d["eft_a_grid"]
    br_ratio = d["eft_grhoa2"] / (d["eft_grhom_t"] * ea**2)
    sel = [0, len(a) // 2, len(a) - 1]
    print(f"{tag} {st['status']}: derivs grhoa2/(grhom_t a^2) at a={a[sel]} -> {ev_ratio[sel]} "
          f"(a^-6 = {a[sel]**-6.0})")
    sel = [0, 200, 400]
    print(f"{tag}: bridge grhoa2/(grhom_t a^2) at a={ea[sel]} -> {br_ratio[sel]}")

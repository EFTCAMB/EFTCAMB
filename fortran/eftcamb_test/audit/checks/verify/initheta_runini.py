"""Verify finding ini-theta-h0-solved-without-eft on the CAMB_RunFromIni path (camb.run_ini == Fortran camb executable).

Writes an ini with thetastar, runs camb.run_ini in a subprocess (A0 tree) and compares the TT output with runs at
hubble = correct H0 (58.86797, EFT-aware brentq) and hubble = the LCDM-solved H0 (63.82937).
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path("/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB")
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")
OUT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit/work/verify/initheta/runini")
OUT.mkdir(parents=True, exist_ok=True)
PARAMS = REPO / "fortran/eftcamb_test/parameters"
THETA = 0.0104110

CHILD = f"""
import sys
sys.path.insert(0, {TREE!r})
import camb
assert camb.__file__.startswith({TREE!r})
camb.run_ini(sys.argv[1])
"""

FLUID = {"EFTflag": "0", "dark_energy_model": "fluid", "w": "-0.8"}
RUNS = [("2_PEFT_wDE_wCDM_1", "thetastar", {}), ("2_PEFT_wDE_wCDM_1", "hubble=58.86797", {}),
        ("2_PEFT_wDE_wCDM_1", "hubble=63.82937", {}), ("1_EFT_GR", "thetastar", FLUID),
        ("1_EFT_GR", "hubble=58.86798", FLUID)]
for model, key, extra in RUNS:
    keys = {"get_scalar_cls": "T", "get_tensor_cls": "F", "get_transfer": "F", "do_lensing": "F",
            "l_max_scalar": "600", "feedback_level": "1", "EFTCAMB_feedback_level": "0",
            "EFTCAMB_write_background": "F", "output_root": str(OUT / f"{model}_{key.replace('=', '')}")} | extra
    # flatten (ini DEFAULT keeps 'hubble', which would conflict with thetastar), so write the full key set
    flat = harness.parse_ini(PARAMS / f"{model}.ini") | keys
    flat.pop("hubble", None)
    if key == "thetastar":
        flat["thetastar"] = str(THETA)
    else:
        flat["hubble"] = key.split("=")[1]
    ini = OUT / f"{model}_{key.replace('=', '')}_in.ini"
    ini.write_text("\n".join(f"{k} = {v}" for k, v in flat.items()) + "\n")
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-c", CHILD, str(ini)], capture_output=True, text=True, cwd=OUT,
                           env=os.environ | {"OMP_NUM_THREADS": "1"}, timeout=600)
    print(f"run_ini {model:20s} {key:16s} {extra} rc={p.returncode} {p.stderr[-300:]}")

# run_ini does not apply feedback_level (global config), so compare the produced TT spectra instead


def tt(model, key):
    return np.loadtxt(OUT / f"{model}_{key.replace('=', '')}_scalCls.dat")[:, 1]


for model, ref in [("2_PEFT_wDE_wCDM_1", ["hubble=58.86797", "hubble=63.82937"]), ("1_EFT_GR", ["hubble=58.86798"])]:
    a = tt(model, "thetastar")
    for r in ref:
        b = tt(model, r)
        print(f"{model:20s} TT(thetastar run) vs TT({r}): max|rel diff| l=2..600 = {np.max(np.abs(a / b - 1)):.3e}")

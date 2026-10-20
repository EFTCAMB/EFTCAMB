"""Repro (track 10a): EFTCAMB forces IntTolBoost >= 2 for every EFTflag != 0 run (results.f90:394-399).

Consequence: an EFT model that is exactly LCDM (2_PEFT_wDE_LCDM_1: pure EFT, wDE = LCDM, no EFT functions) no
longer reproduces the EFTflag = 0 run at default settings; the difference is entirely the silently changed
integration/recombination tolerances (rerunning GR with IntTolBoost = 2 makes the two agree to ~1e-13).
The user's IntTolBoost (e.g. 1) is also silently overridden.

Run: /Users/marco/Software/miniforge3/bin/python t10a-forced-inttolboost-gr-limit.py [TREE]
TREE defaults to the frozen A0 build.
"""

import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
TREE = sys.argv[1] if len(sys.argv) > 1 else (
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0")
sys.path.insert(0, TREE)
import camb  # noqa: E402
import numpy as np  # noqa: E402

PAR = Path(__file__).resolve().parents[2] / "parameters"


def run(ini, itb=None):
    cp = camb.read_ini(str(PAR / ini), no_validate=True)
    if itb is not None:
        cp.Accuracy.IntTolBoost = itb
    r = camb.get_results(cp)
    cl = r.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)
    return cl["tensor"][2:2501], cl["lensed_scalar"][2:2501]


def rel(a, b, col):
    m = b[:, col] != 0
    return np.max(np.abs(a[m, col] / b[m, col] - 1))


t_eft1, s_eft1 = run("2_PEFT_wDE_LCDM_1.ini", itb=1.0)
t_eft2, s_eft2 = run("2_PEFT_wDE_LCDM_1.ini", itb=2.0)
t_gr1, s_gr1 = run("1_EFT_GR.ini", itb=1.0)
t_gr2, s_gr2 = run("1_EFT_GR.ini", itb=2.0)

print("EFT run requested IntTolBoost=1 vs requested IntTolBoost=2 (identical => the request is overridden):")
print(f"  tensor max|diff| rel = {max(rel(t_eft1, t_eft2, c) for c in (0, 1, 2)):.3e}, "
      f"lensed = {max(rel(s_eft1, s_eft2, c) for c in (0, 1, 2)):.3e}")
print("EFT(LCDM-equivalent, wDE=LCDM, no EFT functions) vs EFTflag=0 GR at default settings (expected ~1e-10):")
print(f"  tensor TT max rel diff = {rel(t_eft1, t_gr1, 0):.3e}   tensor BB = {rel(t_eft1, t_gr1, 2):.3e}")
print(f"  lensed TT max rel diff = {rel(s_eft1, s_gr1, 0):.3e}   lensed EE = {rel(s_eft1, s_gr1, 1):.3e}")
print("Same EFT run vs EFTflag=0 GR run with IntTolBoost=2:")
print(f"  tensor TT max rel diff = {rel(t_eft1, t_gr2, 0):.3e}   tensor BB = {rel(t_eft1, t_gr2, 2):.3e}")
print(f"  lensed TT max rel diff = {rel(s_eft1, s_gr2, 0):.3e}   lensed EE = {rel(s_eft1, s_gr2, 1):.3e}")

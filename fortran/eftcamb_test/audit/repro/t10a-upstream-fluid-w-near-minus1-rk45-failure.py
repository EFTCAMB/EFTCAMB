"""Repro (track 10a, upstream CAMB v2): CAMB's own fluid dark energy (EFTflag = 0, cs2 = 1) with a tabulated w(a)
that stays within ~2e-7 of w = -1 during matter domination makes the new RungeKuttaDP45 scalar integrator abort
("RungeKuttaDP45 error -3: ... step-size that is less than or equal to hmin") in plain CAMB v2 (B0) and in the code
under audit (A0, EFTflag = 0).  CAMB 1.6.5 + EFTCAMB (B1, dverk integrator) runs the same input.

The w(a) table is the background equation of state of the 5_quint_cubic test model (A0, eft_gpresq/eft_grhoq on
a = logspace(-8, 0, 401)), stored in fortran/eftcamb_test/audit/inis/t10a_w_5_quint_cubic.dat; the input is otherwise
the EFTCAMB test base_params.ini (CMB, lensing, two number-count/lensing source windows, P(k)).  Clipping the table to 1+w >= 1e-5 (t10a_w_5_quint_cubic_clip.dat) makes A0/B0 run.

Expected: every tree computes the spectra.  Observed (2026-10-17): B0 and A0 fail in calc_transfer with IntTolBoost
1 and 2; B1 runs with IntTolBoost = 1 but fails the same way with IntTolBoost = 2; the clipped table runs everywhere.
The fluid velocity equation is singular as 1+w -> 0 (CAMB recommends PPF there), so this is an upstream robustness
limit, not an EFTCAMB defect.

Run: /Users/marco/Software/miniforge3/bin/python t10a-upstream-fluid-w-near-minus1-rk45-failure.py [TREE ...]
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
BASE = AUDIT.parent / "parameters" / "base_params.ini"
ROOT = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/")
TREES = sys.argv[1:] or [ROOT + t for t in ("B0", "A0", "B1")]

CODE = r"""
import sys; sys.path.insert(0, sys.argv[1]); import camb, numpy as np
cp = camb.read_ini(sys.argv[2], no_validate=True)
cp.Accuracy.IntTolBoost = float(sys.argv[3])
r = camb.get_results(cp)
cl = r.get_cmb_power_spectra(CMB_unit="muK")["total"]
print("ok  D_TT(220) = %.3f" % cl[220, 0])
"""

with tempfile.TemporaryDirectory() as tmp:
    for table in ("t10a_w_5_quint_cubic.dat", "t10a_w_5_quint_cubic_clip.dat"):
        wtab = AUDIT / "inis" / table
        ini = Path(tmp) / f"{table}.ini"
        ini.write_text(f"EFTflag = 0\nuse_tabulated_w = T\nwafile = {wtab}\ncs2_lam = 1\n"
                       f"output_root = {tmp}/x\nEFTCAMB_write_background = F\nDEFAULT({BASE})\n")
        import numpy as np

        w = np.loadtxt(wtab)[:, 1]
        print(f"table {table}: min(1+w) = {np.min(1 + w):.1e}")
        for tree in TREES:
            for itb in ("1", "2"):
                p = subprocess.run([sys.executable, "-c", CODE, tree, str(ini), itb], capture_output=True, text=True,
                                   env=dict(os.environ, OMP_NUM_THREADS="2"), cwd=tmp)
                out = (p.stdout.strip().splitlines() or [""])[-1]
                err = [x.strip() for x in (p.stderr + p.stdout).splitlines() if "RungeKutta" in x or "rror" in x][:1]
                status = out if p.returncode == 0 else f"FAIL rc={p.returncode} {err[0][:80] if err else ''}"
                print(f"  {Path(tree).name:3s} IntTolBoost={itb} -> {status}")

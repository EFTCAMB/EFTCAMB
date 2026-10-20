"""Repro (track 10a, upstream CAMB): camb.get_background(params, no_thermo=False) crashes with SIGSEGV when the
parameters carry source windows (num_redshiftwindows > 0, as in every EFTCAMB test ini through base_params.ini) and
WantCls = False.  Same crash in plain CAMB v2 (B0), pre-rebase EFTCAMB (B1) and the code under audit (A0).

Run: /Users/marco/Software/miniforge3/bin/python t10a-upstream-get-background-windows-segfault.py [TREE ...]
Each case runs in a subprocess; the expected result is "ok" for all of them.
"""

import os
import subprocess
import sys
from pathlib import Path

ROOT = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/")
TREES = sys.argv[1:] or [ROOT + t for t in ("A0", "B0", "B1")]
INI = Path(__file__).resolve().parents[2] / "parameters" / "1_EFT_GR.ini"

CODE = r"""
import sys; sys.path.insert(0, sys.argv[1]); import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
case = int(sys.argv[3])
if case >= 1:
    cp.WantCls = False; cp.WantTransfer = False
if case == 2:
    cp.SourceWindows = []
r = camb.get_background(cp, no_thermo=False)
print('zstar', r.get_derived_params()['zstar'])
"""

labels = {0: "WantCls=T (as read)", 1: "WantCls=F, windows kept", 2: "WantCls=F, windows removed"}
for tree in TREES:
    for case in (0, 1, 2):
        p = subprocess.run([sys.executable, "-c", CODE, tree, str(INI), str(case)], capture_output=True, text=True,
                           env=dict(os.environ, OMP_NUM_THREADS="1"))
        status = "ok" if p.returncode == 0 else f"CRASH rc={p.returncode}"
        print(f"{Path(tree).name:4s} {labels[case]:28s} -> {status}")

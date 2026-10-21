"""Verifier: does upstream get_time_evolution after transfer-only run crash/garbage in B0 and A0 (GR, no EFT)?"""
import subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness
W = r"""
import os, sys
os.environ["OMP_NUM_THREADS"]="1"
sys.path.insert(0, sys.argv[1])
import camb, numpy as np
cp = camb.set_params(H0=67, ombh2=0.022, omch2=0.12, As=2e-9, ns=0.96)
cp.WantCls = sys.argv[2]=="cls"; cp.WantScalars = cp.WantCls; cp.WantTransfer=True; cp.max_l=300
r = camb.get_results(cp)
print("WantScalars in results:", r.Params.WantScalars, flush=True)
ev = r.get_time_evolution([0.05], np.array([100., 1000., 10000.]), ["delta_cdm","Weyl"])
print("ok", ev[0].tolist(), flush=True)
"""
for tree in ["B0", "A0"]:
    for mode in ["cls", "transfer"]:
        with harness.slot(1):
            p = subprocess.run([sys.executable, "-I", "-c", W, str(harness.tree_path(tree)), mode], capture_output=True, text=True, timeout=600)
        print(tree, mode, "rc=", p.returncode, p.stdout.strip().replace("\n", " | "), p.stderr.strip()[-200:])

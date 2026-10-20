"""Run the t9 repro pert-turn-on-nonpositive-silently-gr.py unchanged, under one machine-wide CPU slot at nice 10.

Usage: python pertturnon_run_repro.py [TREE] > out.txt
"""

import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REPRO = HERE.parents[1] / "repro" / "pert-turn-on-nonpositive-silently-gr.py"
args = [sys.executable, str(REPRO)] + sys.argv[1:]
env = dict(os.environ, OMP_NUM_THREADS="1")
with harness.slot(1):
    r = subprocess.run(["nice", "-n", "10"] + args, env=env, timeout=1500)
sys.exit(r.returncode)

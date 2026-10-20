"""Run a script (argv[1:]) as a subprocess while holding harness CPU slots, at nice 10, OMP_NUM_THREADS=1.

Usage: python hderivs_run_slot.py <nslots> <script> [args...]
"""

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

n = int(sys.argv[1])
env = dict(os.environ, OMP_NUM_THREADS=str(n))
env.pop("PYTHONPATH", None)
with harness.slot(n):
    p = subprocess.run([sys.executable, *sys.argv[2:]], env=env, preexec_fn=lambda: os.nice(10))
sys.exit(p.returncode)

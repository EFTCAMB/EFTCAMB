"""Run a command (e.g. a standalone repro script) at nice 10 while holding N machine-wide CPU slots.

usage: python run_with_slot.py N cmd [args...]
"""

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness  # noqa: E402

n = int(sys.argv[1])
env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", str(n)))
with harness.slot(n):
    rc = subprocess.run(sys.argv[2:], env=env, preexec_fn=lambda: os.nice(10)).returncode
sys.exit(rc)

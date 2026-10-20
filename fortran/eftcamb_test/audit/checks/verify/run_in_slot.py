"""Run a python script (args) holding one harness CPU slot, OMP_NUM_THREADS=1, nice 10."""
import os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import harness
env = dict(os.environ, OMP_NUM_THREADS="1")
with harness.slot(1):
    sys.exit(subprocess.run(["nice", "-n", "10", sys.executable] + sys.argv[1:], env=env).returncode)

"""Run a python script (argv[1:]) as a subprocess while holding one harness CPU slot."""
import os, subprocess, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import harness
with harness.slot(1):
    sys.exit(subprocess.call([sys.executable] + sys.argv[1:], env=dict(os.environ, OMP_NUM_THREADS="1")))

"""Run a command while holding n harness CPU slots: run_slotted.py N cmd..."""
import subprocess, sys
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
n = int(sys.argv[1])
with harness.slot(n):
    sys.exit(subprocess.run(["nice", "-n", "10"] + sys.argv[2:]).returncode)

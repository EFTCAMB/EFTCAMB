import sys, runpy
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
with harness.slot(1):
    sys.argv = [sys.argv[1]]
    runpy.run_path(sys.argv[0], run_name="__main__")

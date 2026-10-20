"""Run repro fr-rgr-threshold-high-k-pk with repo path fixed (repo moved from ~/Desktop/IW to ~/IW), inside a CPU slot."""
import importlib.util
import sys

R = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit"
sys.path.insert(0, R + "/checks")
import harness  # noqa: E402

spec = importlib.util.spec_from_file_location("repro", R + "/repro/fr-rgr-threshold-high-k-pk.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
m.PARAMS = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"
with harness.slot(2):
    m.main()

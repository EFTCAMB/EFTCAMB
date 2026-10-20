"""Run the fr-rgr-threshold-high-k-pk repro unchanged except for PARAMS (repo moved from ~/Desktop/IW to ~/IW),
holding a harness CPU slot for its 2-thread child runs."""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

spec = importlib.util.spec_from_file_location("repro", HERE.parents[1] / "repro" / "fr-rgr-threshold-high-k-pk.py")
repro = importlib.util.module_from_spec(spec)
spec.loader.exec_module(repro)
repro.PARAMS = str(harness.PARAMS_DIR)
with harness.slot(2):
    repro.main()

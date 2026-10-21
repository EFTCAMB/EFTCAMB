"""Run the t3-eft-curvature-ignored repro unchanged, but write outputs to the verify work dir (do not clobber t3-gr)."""

import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
spec = importlib.util.spec_from_file_location(
    "repro", REPO / "fortran" / "eftcamb_test" / "audit" / "repro" / "t3-eft-curvature-ignored.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
mod.OUT = mod.harness.AUDIT_ROOT / "work" / "verify" / "t3curv_repro"
mod.OUT.mkdir(parents=True, exist_ok=True)
print("REPO", REPO, "OUT", mod.OUT, file=sys.stderr)
mod.main()

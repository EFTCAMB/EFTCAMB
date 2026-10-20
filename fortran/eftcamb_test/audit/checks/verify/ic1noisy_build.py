"""Verifier (ic1-pidot-transient-noisy-coefficients): build the instrumented A0 copy (IC dump) with release flags."""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

dest = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy" / "A0_icdump"
print(harness.build_tree("A0", dest))

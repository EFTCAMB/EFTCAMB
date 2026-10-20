"""Run a repro script under a harness CPU slot at nice 10 (verify track helper)."""

import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import harness  # noqa: E402

if __name__ == "__main__":
    n = int(os.environ.get("SLOT_N", "2"))
    with harness.slot(n):
        r = subprocess.run(["nice", "-n", "10", sys.executable] + sys.argv[1:])
    sys.exit(r.returncode)

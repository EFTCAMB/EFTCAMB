"""Run the ic23-stale-cache-gammas repro script exactly as given, inside a harness CPU slot (1 thread), at nice 10."""

import os
import subprocess
import sys

sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

REPRO = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/repro/ic23-stale-cache-gammas.py"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "staleg" / f"repro_run{sys.argv[1] if len(sys.argv) > 1 else ''}.txt"

with harness.slot(1):
    p = subprocess.run([sys.executable, REPRO], capture_output=True, text=True, preexec_fn=lambda: os.nice(10),
                       cwd="/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB")
OUT.write_text(p.stdout + "\n--- stderr tail ---\n" + p.stderr[-3000:])
print(p.stdout)
print("rc", p.returncode, "stderr tail:", p.stderr[-1500:])

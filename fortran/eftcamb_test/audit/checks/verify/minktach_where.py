"""Where does EFT_minkowski_limit fail for 5_quint_exponential at default EFTCAMB_stability_time (A0 / given tree)?
Usage: python minktach_where.py [tree]"""
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREE = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/trees/A0"
PAR = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"
W = r'''
import sys
sys.path.insert(0, sys.argv[1])
import camb
for m in ("5_quint_exponential", "5_quint_inverse_power"):
    for st in (None, 1e-6, 1e-5, 1e-4, 1e-3):
        cp = camb.read_ini(f"%s/{m}.ini" % sys.argv[2], no_validate=True)
        cp.EFTCAMB.EFT_minkowski_limit = True
        cp.EFTCAMB.EFTCAMB_feedback_level = 3
        cp.EFTCAMB.EFTCAMB_write_background = False
        if st is not None:
            cp.EFTCAMB.EFTCAMB_stability_time = st
        print("CASE", m, "stability_time", st if st else cp.EFTCAMB.EFTCAMB_stability_time, flush=True)
        try:
            camb.get_background(cp); print("VERDICT stable", flush=True)
        except Exception:
            print("VERDICT unstable", flush=True)
'''
env = dict(os.environ, OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
with harness.slot(1):
    p = subprocess.run([sys.executable, "-c", W, TREE, PAR], capture_output=True, text=True, env=env,
                       cwd=f"{ROOT}/work/verify/minktach", preexec_fn=lambda: os.nice(10), timeout=900)
for line in p.stdout.splitlines():
    if line.startswith(("CASE", "VERDICT")) or "Instability detected" in line:
        print(line)

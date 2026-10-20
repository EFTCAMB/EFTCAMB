"""lldb backtrace of the F-020 crashes (one CAMB slot).

Usage: python verify_getters_segv_lldb.py TREE default_num_params|failmodel_get_results|failmodel_get_background
Runs verify_getters_segv.py's CHILD code under lldb in batch mode and prints the crashing frames.
"""

import os
import subprocess
import sys

sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness  # noqa: E402
import verify_getters_segv as v  # noqa: E402

tree, case = sys.argv[1], sys.argv[2]
work = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/" \
       "scratchpad/audit/work/verify/getters_segv"
child = os.path.join(work, "lldb_child.py")
with open(child, "w") as f:
    f.write(v.CHILD)
env = dict(os.environ, OMP_NUM_THREADS="1", VERIFY_CWD=work)
env.pop("PYTHONPATH", None)
cmd = ["lldb", "--batch", "-o", "run", "-o", "bt 8", "-o", "frame select 0", "-o", "register read pc",
       "-k", "bt 8", "--", sys.executable, child, tree, case, v.INI]
with harness.slot(1):
    r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=600)
lines = r.stdout.splitlines()
keep = [ln for ln in lines if ("frame #" in ln or "stop reason" in ln or "STATE" in ln or "VALUE" in ln
                              or "pc =" in ln or "initialize_parameters" in ln)]
print("\n".join(ln[:220] for ln in keep[:40]))

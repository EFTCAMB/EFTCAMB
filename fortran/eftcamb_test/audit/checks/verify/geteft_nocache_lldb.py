"""Run geteft_nocache_lldb_child.py under lldb (batch) to locate the faulting source line in EFTCAMB_GetEFTFunctions.

argv: tree case [tag]; output to work/verify/geteft_noc/lldb/<tag>_<case>.txt
"""

import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import harness  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/geteft_noc/lldb"
tree, case = sys.argv[1], sys.argv[2]
tag = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(tree)
cmd = ["nice", "-n", "10", "lldb", "--batch", "-o", "run", "-k", "bt 6", "-k", "frame select 0", "-k", "frame variable",
       "-k", "disassemble -p -c 4", "-k", "register read x0 x1 x8 pc", "-k", "kill", "--",
       sys.executable, os.path.join(HERE, "geteft_nocache_lldb_child.py"), tree, case]
import signal  # noqa: E402

with harness.slot(1):
    pr = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True,
                          env=dict(os.environ, OMP_NUM_THREADS="1"))
    try:
        so, se = pr.communicate(timeout=150)
    except subprocess.TimeoutExpired:
        os.killpg(pr.pid, signal.SIGKILL)
        so, se = pr.communicate()
        so += "\nTIMEOUT"
    r = subprocess.CompletedProcess(cmd, pr.returncode, so, se)
os.makedirs(OUT, exist_ok=True)
path = os.path.join(OUT, f"{tag}_{case}.txt")
with open(path, "w") as f:
    f.write(r.stdout + "\n--- stderr ---\n" + r.stderr)
keep = [ln for ln in r.stdout.splitlines() if any(s in ln for s in ("CHILD", "stop reason", "frame #", "EXC_", "at 11_",
                                                                     "->", "x0 =", "x1 =", "x8 =", "pc ="))]
print("\n".join(keep[:40]))
print("full output:", path)

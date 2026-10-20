"""Release-build probe for finding uninit-cache-fields-garbage: what does the stability NaN check actually read in
EFTcdotdot / EFTGamma1PP? Runs uninit_lldb_child.py (sequence of models in one process) under lldb with
uninit_lldb_probe.py, holding one harness CPU slot.

Usage: python verify_uninit_lldb.py [tree] [model ...]
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
DEFAULT_MODELS = ["5_quint_exponential", "4_Des_F_R_LCDM_1", "5_Kmouflage_3", "5_ADE_1", "5_scg", "5_Horava_eta_1",
                  "5_Kmimic_1", "4_Des_mc5e_wCDM_1", "5_quint_exponential", "5_quint_axion"]

arg = sys.argv[1] if len(sys.argv) > 1 else "A0"
tree = arg if os.path.isabs(arg) else str(ROOT / "trees" / arg)
models = sys.argv[2:] or DEFAULT_MODELS
env = dict(os.environ, OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
with tempfile.TemporaryDirectory() as d, harness.slot(1):
    cmd = ["nice", "-n", "10", "lldb", "--batch", "-o", f"command script import {HERE / 'uninit_lldb_probe.py'}",
           "-o", "run", "-k", "bt 8", "--", sys.executable, str(HERE / "uninit_lldb_child.py"), tree, d, *models]
    p = subprocess.run(cmd, env=env, cwd=d, capture_output=True, text=True, timeout=580)
for ln in p.stdout.splitlines():
    if ln.startswith(("PROBE", "MODEL", "  accepted", "  ERROR", "Process", "error")):
        print(ln)
print("rc", p.returncode, p.stderr[-500:])

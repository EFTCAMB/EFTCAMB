"""Backtrace + local variables at the AD bounds-check abort for 3_RPH_KT_power_1 (lldb, breakpoint on gfortran's
runtime-error entry). Prints the call chain and the x/ind values in the precompute frame."""

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness  # noqa: E402

AUDIT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREE = f"{AUDIT}/trees/AD"
INI = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/" + (sys.argv[1] if len(sys.argv) > 1 else "3_RPH_KT_power_1.ini")
OUT = Path(f"{AUDIT}/work/verify/interp_xfinal/lldb")
child = OUT / "child.py"
child.write_text(
    "import sys\nsys.path.insert(0, sys.argv[1])\nimport camb\n"
    "cp = camb.read_ini(sys.argv[2], no_validate=True)\ncp.WantCls = False\ncp.WantTransfer = False\n"
    "camb.get_background(cp)\nprint('passed')\n"
)
cmds = [
    "settings set target.process.stop-on-exec false",
    "breakpoint set -r _gfortran_runtime_error_at",
    "process launch",
    "bt 12",
    "frame select 1",
    "frame variable x ind",
    "expression -- self->_data->x_initial",
    "expression -- self->_data->x_final",
    "expression -- self->_data->grid_width",
    "expression -- self->_data->num_points",
    "frame select 3",
    "frame variable",
    "kill",
    "quit",
]
args = ["lldb", "-b"]
for c in cmds:
    args += ["-o", c]
args += ["--", sys.executable, str(child), TREE, INI]
env = dict(os.environ, OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
with harness.slot(1):
    r = subprocess.run(args, capture_output=True, text=True, cwd=OUT, env=env, timeout=500)
txt = r.stdout + r.stderr
(OUT / "lldb_out.txt").write_text(txt)
keep = [ln for ln in txt.splitlines() if ln.strip().startswith(("frame #", "*", "(", "x", "ind")) or "= " in ln]
print("\n".join(keep[:80]))

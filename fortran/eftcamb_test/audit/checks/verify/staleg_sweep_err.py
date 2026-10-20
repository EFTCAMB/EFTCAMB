"""Rerun failing sweep entries (staleg_sweep.py) one at a time; show return code / last output lines.
usage: python staleg_sweep_err.py IC model [model ...]"""

import os
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "staleg_sweep.py")).read()
CHILD = re.search(r'CHILD = r"""(.*?)"""', src, re.S).group(1)
TREE = str(harness.TREES["A0"])
OUT = harness.AUDIT_ROOT / "work" / "verify" / "staleg" / "sweep"
ic = sys.argv[1]
for model in sys.argv[2:]:
    tmp = tempfile.mkdtemp(dir=OUT)
    ini = os.path.join(tmp, "w.ini")
    with open(ini, "w") as f:
        f.write(f"DEFAULT({harness.PARAMS_DIR}/{model}.ini)\noutput_root = {tmp}/r\nfeedback_level = 0\n"
                f"EFTCAMB_write_background = F\nEFT_IC_type = {ic}\n")
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-c", CHILD, TREE, ini, "0.05"], capture_output=True, text=True, cwd=tmp,
                           env=env, preexec_fn=lambda: os.nice(10), timeout=600)
    lines = [x for x in (p.stdout + p.stderr).splitlines() if x.strip()]
    print(f"{model} IC{ic} rc={p.returncode}: " + (" | ".join(x[:110] for x in lines[-3:]) if lines else "<no output>"))

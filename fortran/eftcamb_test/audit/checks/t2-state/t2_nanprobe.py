"""Run get_background (SetParams: designer/full background solve, RGR, stability check) for a list of models on a
tree, one fresh process per model, and report status + which timestep-cache field the stability check found NaN
(EFTCAMB_feedback_level = 3 prints '<name> is NaN').  Used with the NaN-instrumented tree (work/t2-state/nanit:
-finit-real=snan -finit-derived + interpolation tables poisoned with NaN at allocation) vs A0."""

from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

CHILD = r'''
import sys
sys.path.insert(0, sys.argv[1])
import camb
open("w.ini", "w").write("get_transfer = F\nget_tensor_cls = F\ndo_lensing = F\nl_max_scalar = 100\n"
                         "feedback_level = 1\nEFTCAMB_feedback_level = 3\n" + "".join(a + "\n" for a in sys.argv[3:])
                         + "DEFAULT(" + sys.argv[2] + ")\n")
cp = camb.read_ini("w.ini")
try:
    r = camb.get_background(cp)
    print("RESULT ok H(z=1)=%r" % r.hubble_parameter(1.0))
except Exception as e:
    print("RESULT error " + str(e).replace("\n", " "))
'''


def one(tree: str, model: str, extra: list[str]) -> str:
    ini = harness.PARAMS_DIR / f"{model}.ini"
    if not ini.exists():
        ini = HERE.parents[1] / "inis" / f"{model}.ini"
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1), tempfile.TemporaryDirectory() as d:
        p = subprocess.run([sys.executable, "-c", CHILD, str(harness.tree_path(tree)), str(ini), *extra], cwd=d,
                           env=env, capture_output=True, text=True, timeout=540, preexec_fn=harness._nice)
    out = p.stdout + p.stderr
    res = [ln for ln in out.splitlines() if ln.startswith("RESULT")]
    nans = sorted(set(re.findall(r"^\s*(\S+)\s+is NaN", out, re.M)))
    inst = re.findall(r"Instability detected at a =\s*(\S+)", out)
    msg = res[-1][7:] if res else f"rc={p.returncode} " + out.strip().splitlines()[-1][:150] if out.strip() else ""
    return f"{model:28s} {msg[:90]:90s} NaN fields: {','.join(nans) or '-'}  unstable at a={inst[:1]}"


if __name__ == "__main__":
    tree = sys.argv[1]
    models = sys.argv[2].split(",")
    extra = sys.argv[3:]
    with ThreadPoolExecutor(int(os.environ.get("JOBS", "2"))) as ex:
        for line in ex.map(lambda m: one(tree, m, extra), models):
            print(line, flush=True)

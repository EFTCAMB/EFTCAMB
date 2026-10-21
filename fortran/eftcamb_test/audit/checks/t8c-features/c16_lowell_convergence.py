"""Check 16 (track t8c-features): low-ell (ISW) TT convergence of 2_PEFT_gamma3_const_1 vs GR, l_max_scalar = 300.

Finding of this check (not filed as a defect): the gamma3 model's unlensed TT at l <~ 40 is not monotonic in the time
sampling: default, TimeStepBoost = 2 / accuracy_boost = 2 and accuracy_boost = 3 differ by up to 1% (accuracy_boost=2
is the outlier), while GR converges monotonically at l_max_scalar = 2500 (check c14 extension).  IntTolBoost = 4
changes nothing.  Same magnitude in the pre-rebase EFTCAMB (B1, check c14/c15), so it is not a CAMB v2 accuracy
default; it matches the source-time-grid dependence of the pi field reported by track t10b-convergence.
Run: /Users/marco/Software/miniforge3/bin/python c16_lowell_convergence.py   (10 short CAMB runs)
"""

import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
try:  # hold machine-wide CPU slots of the audit harness when available
    from harness import slot
except Exception:  # noqa: BLE001
    from contextlib import nullcontext

    def slot(n):
        return nullcontext()


TREE = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
COMMON = ("get_tensor_cls = F\nget_transfer = F\nnum_redshiftwindows = 0\nfeedback_level = 0\n"
          "EFTCAMB_write_background = F\nl_max_scalar = 300\ndo_lensing = F\n")
WORKER = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2])
cl = camb.get_results(cp).get_unlensed_scalar_cls(CMB_unit="muK", raw_cl=True, lmax=40)[:, 0]
print(json.dumps(cl.tolist()))
"""
VARIANTS = {"default": "", "TimeStepBoost=2": "TimeStepBoost = 2\n", "IntTolBoost=4": "IntTolBoost = 4\n",
            "accuracy_boost=2": "accuracy_boost = 2\n", "accuracy_boost=3": "accuracy_boost = 3\n"}
tmp = Path(tempfile.mkdtemp())


def run(job):
    model, tag = job
    ini = tmp / f"{model}_{tag.replace('=', '')}.ini"
    ini.write_text(f"DEFAULT({PARAMS / model}.ini)\n{COMMON}{VARIANTS[tag]}")
    with slot(1):
        p = subprocess.run([sys.executable, "-c", WORKER, TREE, str(ini)], capture_output=True, text=True, cwd=tmp,
                           env=dict(os.environ, OMP_NUM_THREADS="1"), preexec_fn=lambda: os.nice(10))
    lines = [x for x in p.stdout.splitlines() if x.startswith("[")]
    if not lines:
        raise RuntimeError(p.stderr[-800:])
    return np.array(json.loads(lines[-1]))


jobs = [(m, t) for m in ("2_PEFT_gamma3_const_1", "1_EFT_GR") for t in VARIANTS]
with ThreadPoolExecutor(4) as ex:
    cls = dict(zip(jobs, ex.map(run, jobs)))
for m in ("2_PEFT_gamma3_const_1", "1_EFT_GR"):
    ref = cls[(m, "accuracy_boost=3")][2:]
    print(f"{m}: unlensed TT / TT(accuracy_boost=3) - 1, l = 2..40")
    for t in VARIANTS:
        r = cls[(m, t)][2:] / ref - 1
        i = int(np.argmax(np.abs(r)))
        print(f"   {t:18s} max |.| = {abs(r[i]):.2e} at l={i + 2:2d}   l=2: {r[0]:+.2e}  l=13: {r[11]:+.2e}  "
              f"l=30: {r[28]:+.2e}")

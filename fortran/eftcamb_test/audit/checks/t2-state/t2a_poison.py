"""Track 2a: localise the thread-count dependence of get_scale_evolution (008p5 Beyond Horndeski).

Runs get_background + get_scale_evolution in a child process on an instrumented tree (work/t2-state/fix or fixdbg,
see 11_EFTCAMB_python.f90 there): EFT_AUDIT_EVINIT=1 zero-initialises the OpenMP-private EV%eft_cache before
each k, EFT_AUDIT_EVINIT=2 fills it with signalling NaN (with -ffpe-trap=invalid the first arithmetic read of a
field that was never written in this k traps and gfortran prints a backtrace).

Usage: python t2a_poison.py <tree> <model> <evinit 0|1|2> [k1,k2,...]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

CHILD = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2])
r = camb.get_background(cp)
k = [float(x) for x in sys.argv[3].split(",")]
a = np.logspace(-6, 0, 300)
fields, ev = cp.EFTCAMB.get_scale_evolution(r, k, a)
out = {f: np.asarray(ev[f]).tolist() for f in ev.dtype.names}
print("RESULT " + json.dumps(out))
"""


def ini_for(model: str) -> Path:
    p = harness.PARAMS_DIR / f"{model}.ini"
    return p if p.exists() else HERE.parents[1] / "inis" / f"{model}.ini"


def run(tree: str, model: str, evinit: str, ks: str, threads: int = 1) -> tuple[int, str]:
    with harness.slot(threads), tempfile.TemporaryDirectory() as d:
        w = Path(d) / "w.ini"
        w.write_text(f"DEFAULT({ini_for(model)})\noutput_root = {d}/r\nEFTCAMB_write_background = F\n"
                     "feedback_level = 0\nEFTCAMB_feedback_level = 0\n")
        env = dict(os.environ, OMP_NUM_THREADS=str(threads), EFT_AUDIT_EVINIT=evinit)
        env.pop("PYTHONPATH", None)
        p = subprocess.run([sys.executable, "-c", CHILD, str(harness.tree_path(tree)), str(w), ks], cwd=d, env=env,
                           capture_output=True, text=True, timeout=540, preexec_fn=harness._nice)
        return p.returncode, p.stdout + p.stderr


if __name__ == "__main__":
    tree, model, evinit = sys.argv[1], sys.argv[2], sys.argv[3]
    ks = sys.argv[4] if len(sys.argv) > 4 else "1e-3"
    threads = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    rc, out = run(tree, model, evinit, ks, threads)
    lines = [ln for ln in out.splitlines() if not ln.startswith("RESULT ")]
    res = [ln for ln in out.splitlines() if ln.startswith("RESULT ")]
    print(f"rc={rc}")
    print("\n".join(ln for ln in lines if ln.strip() and "1.00000" not in ln[:20])[-6000:])
    if res:
        d = json.loads(res[-1][7:])
        Path(os.environ.get("T2A_POISON_OUT", "/dev/null")).write_text(json.dumps(d)) if os.environ.get(
            "T2A_POISON_OUT") else None
        print("RESULT clxc[k0][:3] =", d["clxc"][0][:3])

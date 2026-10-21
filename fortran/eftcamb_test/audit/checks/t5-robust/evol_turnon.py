"""Track 5: pi(a), pi'(a) and the pi-equation coefficients right after the perturbation turn-on, for a few k.

Runs one model on a tree (P(k)-only settings unless --full), then calls EFTCAMB.get_scale_evolution on a
log-spaced grid starting at the turn-on time.  Prints, for each k, the IC denominator C + k^2 D1 + k^4 D2 at
the first sample, pi, pi' and the quasi-static value -h0 E/(C + k^2 D1 + k^4 D2) at a few times.

Usage: python evol_turnon.py MODEL --k 0.0099,0.0104,0.0109 [--tree A0] [--a0 1e-8] [--a1 1e-5] [--set key=val ...]
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

CHILD = r"""
import sys, json, numpy as np
tree, ini, ks, a0, a1 = sys.argv[1], sys.argv[2], [float(x) for x in sys.argv[3].split(',')], float(sys.argv[4]), float(sys.argv[5])
sys.path.insert(0, tree)
import camb
cp = camb.read_ini(ini, no_validate=True)
res = camb.get_results(cp)
a = np.logspace(np.log10(a0), np.log10(a1), 120)
fields, arr = cp.EFTCAMB.get_scale_evolution(res, ks, a)
h0 = cp.H0 / 299792.458
for ik, k in enumerate(ks):
    C, D1, D2, E = (arr['EFTpi' + x][ik] for x in ('C', 'D1', 'D2', 'E'))
    den = C + k**2 * D1 + k**4 * D2
    print('k=%.6g  first sample a=%.3e: C=%.3e D1=%.3e D2=%.3e den=%.3e E=%.3e' % (k, arr['a'][ik][0], C[0], D1[0], D2[0], den[0], E[0]))
    for i in list(range(0, 10)) + list(range(10, 120, 12)):
        qs = -h0 * E[i] / den[i] if den[i] != 0 else float('nan')
        print('   a=%.3e pi=% .4e pidot=% .4e piQS=% .4e den=% .3e' % (arr['a'][ik][i], arr['pi'][ik][i], arr['pidot'][ik][i], qs, den[i]))
"""


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("model")
    p.add_argument("--k", required=True)
    p.add_argument("--tree", default="A0")
    p.add_argument("--a0", type=float, default=1e-8)
    p.add_argument("--a1", type=float, default=1e-5)
    p.add_argument("--set", nargs="*", default=[])
    a = p.parse_args()
    tree = harness.AUDIT_ROOT / "trees" / a.tree if "/" not in a.tree else Path(a.tree)
    ini = Path(a.model) if "/" in a.model else harness.PARAMS_DIR / f"{a.model}.ini"
    tmp = Path(tempfile.mkdtemp(dir=harness.AUDIT_ROOT / "work" / "t5-robust"))
    lines = [f"DEFAULT({ini})", f"output_root = {tmp / 'x'}", "EFTCAMB_write_background = F", "feedback_level = 0",
             "get_tensor_cls = F"] + a.set
    (tmp / "w.ini").write_text("\n".join(lines) + "\n")
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    with harness.slot(2):
        r = subprocess.run(["nice", "-n", "10", sys.executable, "-c", CHILD, str(tree), str(tmp / "w.ini"), a.k,
                            str(a.a0), str(a.a1)], capture_output=True, text=True, cwd=tmp, env=env, timeout=540)
    print(r.stdout)
    if r.returncode:
        print(r.stderr[-2000:])


if __name__ == "__main__":
    main()

"""Repro (t9-wrapper): at k ~ 1/Mpc the pi field (and pidot, Psi, mu, ... built from it) returned by
EFTCAMB.get_evolution / get_scale_evolution is dominated by integration error at the default tolerance.

GetEFTOutputEvolutionFork (11_EFTCAMB_python.f90:362) integrates with CAMB's source tolerance
tol1 = base_tol/exp(AccuracyBoost*IntTolBoost - 1). For 2_PEFT_Omega_power_1 at k = 1.317/Mpc, pi(a) at default
tolerance oscillates with amplitude ~2e-2 around the converged solution (~7e-4, obtained with IntTolBoost = 8), and
it changes O(1) when the same times are requested with a different set of output stops (every 4th time). The same
noise is present in the DEBUG source dump of the Cl computation, but the CMB / lensing spectra are insensitive to it
(printed below), so the defect is in the accuracy of what get_evolution returns at high k.
Expected: get_evolution output independent of the requested sampling and converged (|pi_default - pi_tight| small
compared with |pi|).

Run: /Users/marco/Software/miniforge3/bin/python get-evolution-pi-highk-tolerance.py [--tree PATH]
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
INI = REPO / "fortran" / "eftcamb_test" / "parameters" / "2_PEFT_Omega_power_1.ini"
TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")

CHILD = r"""
import sys, os
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
ini, tmp = sys.argv[2], sys.argv[3]
fast = os.path.join(tmp, 'model.ini')
open(fast, 'w').write(f'DEFAULT({ini})\nl_max_scalar = 1500\nget_tensor_cls = F\nget_transfer = F\n'
                      'EFTCAMB_write_background = F\nfeedback_level = 0\n')
k = 1.3172836143
a = np.logspace(np.log10(0.02), 0, 400)
out = {}
for boost in (1.0, 8.0):
    p = camb.read_ini(fast)
    p.Accuracy.IntTolBoost = boost
    r = camb.get_results(p)
    fine = p.EFTCAMB.get_scale_evolution(r, k, a)[1][0]
    coarse = p.EFTCAMB.get_scale_evolution(r, k, a[::4])[1][0]
    pp = r.get_lens_potential_cls(lmax=1500, raw_cl=True)[:, 0]
    out[boost] = (fine, coarse, pp)
f1, c1, pp1 = out[1.0]
f8, c8, pp8 = out[8.0]
print('   a     pi(default)   pi(default, every 4th stop)   pi(IntTolBoost=8)   Psi(default)   Psi(IntTolBoost=8)')
for i in range(0, len(a), 50):
    j = i // 4 * 4
    print(f'{a[j]:6.3f}  {f1["pi"][j]: .4e}     {c1["pi"][j // 4]: .4e}                  {f8["pi"][j]: .4e}'
          f'        {f1["Psi"][j]: .4e}    {f8["Psi"][j]: .4e}')
for name in ('pi', 'pidot', 'Psi', 'mu'):
    m = np.max(np.abs(f8[name]))
    print(f'{name:6s}: max|default - tight|/max|tight| = {np.max(np.abs(f1[name] - f8[name])) / m:9.3e} ;'
          f' default fine-vs-coarse stops = {np.max(np.abs(f1[name][::4] - c1[name])) / m:9.3e} ;'
          f' tight fine-vs-coarse = {np.max(np.abs(f8[name][::4] - c8[name])) / m:9.3e}   (expected << 1)')
print(f'C_L^phiphi (L=2..1500): max|default/tight - 1| = {np.max(np.abs(pp1[2:] / pp8[2:] - 1)):.2e}')
"""

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=TREE)
    a = ap.parse_args()
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run([sys.executable, "-c", CHILD, a.tree, str(INI), tmp], cwd=tmp, env=env,
                           capture_output=True, text=True)
        print("\n".join(ln for ln in r.stdout.splitlines() if not ln.startswith((" *", "*"))))
        if r.returncode:
            print("child failed rc =", r.returncode, r.stderr[-1500:])

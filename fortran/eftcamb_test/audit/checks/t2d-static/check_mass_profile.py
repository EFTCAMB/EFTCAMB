"""Profile of the mass-stability quantities mu_1/H^2, mu_2/H^2 (EFT_mu1*(a/adotoa)^2) versus a, computed by the
public EFTTestStability (so compute_additional_derivs + compute_mass_factors run), on a given tree.

Usage: python check_mass_profile.py TREE MODEL [MODEL ...]   (MODEL = test ini stem)
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
PAR = REPO / "fortran/eftcamb_test/parameters"
CHILD = r"""
import sys, json, ctypes, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
from camb import eftcamb
from camb.baseconfig import camblib
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.EFTCAMB.EFT_mass_stability = False
res = camb.get_background(pars)
P = res.Params
P.EFTCAMB.EFT_mass_stability = True
P.EFTCAMB.EFTCAMB_feedback_level = 0
f = camblib.__eftcamb_stability_MOD_eftteststability
f.restype = ctypes.c_int
ton = P.EFTCAMB.EFTCAMB_pert_turn_on
out = []
for a in np.concatenate([ton * np.array([1.001, 1.01, 1.03, 1.1, 1.3]), np.logspace(np.log10(ton * 2), 0, 12)]):
    c = eftcamb.EFTCAMB_timestep_cache()
    kmax = ctypes.c_double(10.0)
    ok = f(ctypes.byref(ctypes.c_double(a)), ctypes.byref(kmax), ctypes.byref(P.EFTCAMB),
           ctypes.byref(P.EFTCAMB_parameter_cache), ctypes.byref(c))
    h2 = (c.adotoa / a) ** 2
    out.append([float(a), bool(ok), c.EFT_mu1 / h2 * 1.0, c.EFT_mu2 / h2, c.EFTcdotdot, c.EFTLambdadotdot,
                c.EFTGamma1PP, c.EFTGamma2PP])
print("JSON" + json.dumps(dict(ton=ton, rows=out)))
"""


def main():
    tree = Path(sys.argv[1])
    tmp = Path(tempfile.mkdtemp(prefix="mprof_"))
    for m in sys.argv[2:]:
        w = tmp / f"{m}.ini"
        w.write_text(f"DEFAULT({PAR / (m + '.ini')})\nfeedback_level = 0\nEFTCAMB_write_background = F\n"
                     f"output_root = {tmp / m}\n")
        r = subprocess.run([sys.executable, "-c", CHILD, str(tree), str(w)], capture_output=True, text=True, cwd=tmp,
                           env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900)
        js = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
        print(f"== {m} on {tree.name}")
        if not js:
            print(r.stdout[-800:], r.stderr[-800:])
            continue
        d = json.loads(js[-1][4:])
        print(f"   turn-on a = {d['ton']}")
        print("        a     stable    mu1/H^2      mu2/H^2      cdotdot       Lambdadotdot  Gamma1PP     Gamma2PP")
        for a, ok, m1, m2, cdd, ldd, g1, g2 in d["rows"]:
            print(f"   {a:9.3e}  {str(ok):5s} {m1:+.4e} {m2:+.4e} {cdd:+.4e} {ldd:+.4e} {g1:+.4e} {g2:+.4e}")


if __name__ == "__main__":
    main()

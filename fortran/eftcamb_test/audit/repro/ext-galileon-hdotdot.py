"""Repro: Extended Galileon (EFTflag=4, FullMappingEFTmodel=7) returns a wrong d^2H/dtau^2 on the tracker branch.

Expected: the tracker Friedmann equation of arXiv:1912.07586 eq. (16),
    E^(2+s) = Omega_phi0 + Omega_tot(a) E^s,   E = H/(a H0), H = conformal Hubble,
fixes H(a); Hdotdot = d^2H/dtau^2 is obtained here by solving eq. (16) with brentq and differentiating
numerically (d/dtau = H d/dln a).  Third route: finite difference of the code's own Hdot.
Observed: the code's cached Hdotdot (get_eft_functions) is ~1e20 times too large for a > 0.2.

Run: /Users/marco/Software/miniforge3/bin/python ext-galileon-hdotdot.py
"""

import json
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy.optimize import brentq

TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
        "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
BASE = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/base_params.ini"
S, Q = 2.0, 0.35  # GCCG1 of arXiv:1912.07586 Table I

INI = f"""EFTflag = 4
FullMappingEFTmodel = 7
ExtendedGalileon_S = {S}
ExtendedGalileon_q = {Q}
omnuh2 = 0
massive_neutrinos = 0
massless_neutrinos = 3.046
nu_mass_eigenstates = 0
feedback_level = 0
DEFAULT({BASE})
"""

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_background(cp)
a = np.exp(np.linspace(np.log(0.25), 0.0, 41))
_, arr = cp.EFTCAMB.get_eft_functions(res, a)
print('JSON' + json.dumps({'a': a.tolist(), 'H': arr['adotoa'].tolist(), 'Hdot': arr['Hdot'].tolist(),
    'Hdotdot': arr['Hdotdot'].tolist(), 'grhom': (arr['grhob_t'] + arr['grhoc_t']).tolist(),
    'grhor': (arr['grhor_t'] + arr['grhog_t']).tolist(), 'H0': cp.H0}))
"""


def main() -> None:
    fd, ini = tempfile.mkstemp(suffix=".ini")
    os.close(fd)
    with open(ini, "w") as f:
        f.write(INI)
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", WORKER, TREE, ini], capture_output=True, text=True, env=env, timeout=900)
    os.remove(ini)
    js = [ln for ln in p.stdout.splitlines() if ln.startswith("JSON")]
    if not js:
        print(p.stdout[-2000:], p.stderr[-2000:])
        raise SystemExit("CAMB worker failed")
    d = {k: np.array(v) if isinstance(v, list) else v for k, v in json.loads(js[0][4:]).items()}
    a, H0 = d["a"], d["H0"] / 299792.458
    Om = d["grhom"][-1] / (3 * H0**2)  # a = 1 entry
    Orad = d["grhor"][-1] / (3 * H0**2)
    Ophi = 1 - Om - Orad

    def Hconf(x):
        ot = Om / x**3 + Orad / x**4
        e = brentq(lambda y: y ** (2 + S) - Ophi - ot * y**S, 1e-3, 1e6, xtol=1e-15, rtol=1e-15)
        return x * H0 * e

    def dlna(f, x, h=1e-4):
        lx = np.log(x)
        v = [f(np.exp(lx + h * j)) for j in (-2, -1, 1, 2)]
        return (v[0] - 8 * v[1] + 8 * v[2] - v[3]) / (12 * h)

    def Hdot(x):
        return Hconf(x) * dlna(Hconf, x)

    exp_hdd = np.array([Hconf(x) * dlna(Hdot, x) for x in a])
    fd_code = d["H"] * np.gradient(d["Hdot"], np.log(a))  # third route: FD of the code's own Hdot
    print(f"{'a':>6} {'Hdotdot code':>14} {'expected':>14} {'FD(code Hdot)':>14}")
    for i in range(0, len(a), 8):
        print(f"{a[i]:6.3f} {d['Hdotdot'][i]:14.6e} {exp_hdd[i]:14.6e} {fd_code[i]:14.6e}")
    worst = np.max(np.abs(d["Hdotdot"] / exp_hdd - 1))
    worst_fd = np.max(np.abs(fd_code[2:-2] / exp_hdd[2:-2] - 1))
    print(f"observed max |Hdotdot_code/expected - 1| = {worst:.3e}   (expected <~1e-6)")
    print(f"third route max |FD(code Hdot)/expected - 1| = {worst_fd:.3e}  (FD on 41-point grid)")
    print("RESULT:", "BUG REPRODUCED" if worst > 1e-2 else "not reproduced")


if __name__ == "__main__":
    main()

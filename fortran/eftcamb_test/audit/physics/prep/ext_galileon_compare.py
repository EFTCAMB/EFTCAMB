"""Compare EFTCAMB (A0) Extended Galileon background with the independent derivation in ext_galileon_background.py.

Checks (arXiv:1912.07586 eq. 16 tracker Friedmann equation, conformal Hubble H):
  1. residual of eq. 16 for the code's adotoa,
  2. code Hdot    vs implicit-differentiation Hdot  (evaluated at the code's adotoa),
  3. code Hdotdot vs implicit-differentiation Hddot (evaluated at the code's adotoa),
  4. third route: code Hdotdot vs finite difference of the code's own Hdot along conformal time.
Runs CAMB in a subprocess (no massive neutrinos so that Omega_tot = Om/a^3 + Or/a^4 exactly).

Run: /Users/marco/Software/miniforge3/bin/python ext_galileon_compare.py [S q]
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
        "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
INI = HERE.parents[1] / "inis" / "008p7-ext-galileon.ini"

WORKER = r"""
import sys, json, os, tempfile
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
ini, S, q = sys.argv[2], sys.argv[3], sys.argv[4]
fd, tmp = tempfile.mkstemp(suffix='.ini'); os.close(fd)
open(tmp, 'w').write(f'DEFAULT({ini})\nExtendedGalileon_S = {S}\nExtendedGalileon_q = {q}\nomnuh2 = 0\n'
                     'massive_neutrinos = 0\nmassless_neutrinos = 3.046\nnu_mass_eigenstates = 0\nfeedback_level = 0\n')
cp = camb.read_ini(tmp, no_validate=True); os.remove(tmp)
res = camb.get_background(cp)
a = np.concatenate([np.linspace(0.21, 1.0, 80), np.linspace(0.05, 0.199, 30)])
_, arr = cp.EFTCAMB.get_eft_functions(res, a)
pc = {k: getattr(cp, k) for k in ('H0',)}
out = {'a': a.tolist(), 'adotoa': arr['adotoa'].tolist(), 'Hdot': arr['Hdot'].tolist(),
       'Hdotdot': arr['Hdotdot'].tolist(), 'grhob_t': arr['grhob_t'].tolist(), 'grhoc_t': arr['grhoc_t'].tolist(),
       'grhor_t': arr['grhor_t'].tolist(), 'grhog_t': arr['grhog_t'].tolist(), 'H0': cp.H0,
       'omegade': res.get_Omega('de', z=0)}
print('JSON' + json.dumps(out))
"""


def main() -> None:
    S = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
    q = float(sys.argv[2]) if len(sys.argv) > 2 else 0.35
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", WORKER, TREE, str(INI), str(S), str(q)], capture_output=True, text=True,
                       env=env, timeout=900)
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("JSON")]
    if not line:
        print(p.stdout[-3000:], p.stderr[-3000:])
        raise SystemExit("worker failed")
    d = {k: np.array(v) if isinstance(v, list) else v for k, v in json.loads(line[0][4:]).items()}
    sys.path.insert(0, str(HERE))
    from ext_galileon_background import f_Hddot, f_Hdot  # independent expressions

    a, Hc = d["a"], d["adotoa"]
    H0 = d["H0"] / 299792.458  # 1/Mpc
    # Omega_tot(a) from the code's own densities (8 pi G rho a^2 in 1/Mpc^2): Omega_i(a) = grho_i/(3 H0^2 a^2)
    ot_m = (d["grhob_t"] + d["grhoc_t"]) / (3 * H0**2 * a**2)
    ot_r = (d["grhor_t"] + d["grhog_t"]) / (3 * H0**2 * a**2)
    Om, Or = ot_m[0] * a[0] ** 3, ot_r[0] * a[0] ** 4
    Ophi = 1 - Om - Or
    Ot = Om / a**3 + Or / a**4
    Otp = -3 * Om / a**4 - 4 * Or / a**5
    Otpp = 12 * Om / a**5 + 20 * Or / a**6
    E = Hc / (a * H0)
    trk = a > 0.2
    resid = E ** (2 + S) - Ophi - Ot * E**S
    hd_ref = f_Hdot(a, H0, S, Ophi, E, Ot, Otp)
    hdd_ref = f_Hddot(a, H0, S, Ophi, E, Ot, Otp, Otpp)
    r_hd = np.abs(d["Hdot"] / hd_ref - 1)
    r_hdd = np.abs(d["Hdotdot"] / hdd_ref - 1)
    # third route: finite difference of the code's Hdot in conformal time on the tracker branch, d/dtau = H d/dlna
    m = trk
    order = np.argsort(a[m])
    la, hdc, Hcm, hddc = np.log(a[m][order]), d["Hdot"][m][order], Hc[m][order], d["Hdotdot"][m][order]
    fd = Hcm * np.gradient(hdc, la)
    r_fd_code = np.abs(hddc[2:-2] / fd[2:-2] - 1)
    r_fd_ref = np.abs(hdd_ref[m][order][2:-2] / fd[2:-2] - 1)
    print(f"S={S} q={q}  Omega_m0={Om:.5f} Omega_r0={Or:.3e} Omega_phi0={Ophi:.5f}")
    print(f"CHECK 1 [1912.07586 eq.16] tracker residual of code adotoa (a>0.2): max {np.max(np.abs(resid[trk])):.2e}")
    print(f"CHECK 2 [1912.07586 eq.16 d/dtau] code Hdot vs derivation (a>0.2): max rel {np.max(r_hd[trk]):.2e}")
    print(f"CHECK 3 [1912.07586 eq.16 d2/dtau2] code Hdotdot vs derivation (a>0.2): max rel {np.max(r_hdd[trk]):.2e}"
          f"  median {np.median(r_hdd[trk]):.2e}")
    print(f"CHECK 4 [third route] code Hdotdot vs FD of code Hdot: max rel {np.max(r_fd_code):.2e};"
          f" derivation vs FD of code Hdot: max rel {np.max(r_fd_ref):.2e}")
    for i in np.where(trk)[0][::16]:
        print(f"   a={a[i]:.3f} tracker residual={resid[i]: .3e} E={E[i]:.8f}")
    for i in np.where(trk)[0][::16]:
        print(f"   a={a[i]:.3f} Hdotdot code={d['Hdotdot'][i]: .6e} derivation={hdd_ref[i]: .6e}")
    # continuity of adotoa across the a=0.2 switch between Friedmann (a<0.2) and tracker (a>0.2) branches
    lo = np.argmin(np.where(a < 0.2, 0.2 - a, np.inf))
    hi = np.argmin(np.where(a > 0.2, a - 0.2, np.inf))
    Etr_lo = None
    print(f"   adotoa at a={a[lo]:.4f}: {Hc[lo]:.8e}  at a={a[hi]:.4f}: {Hc[hi]:.8e}")
    # tracker solution at the a<0.2 points: solve eq.16 and compare to code
    from scipy.optimize import brentq
    devs = []
    for i in np.where(~trk)[0]:
        e = brentq(lambda x: x ** (2 + S) - Ophi - Ot[i] * x**S, 1e-3, 1e8, xtol=1e-14, rtol=1e-14)
        devs.append(abs(Hc[i] / (a[i] * H0 * e) - 1))
    print(f"CHECK 5 [1912.07586 eq.16] code adotoa vs tracker solution for 0.05<a<0.2: max rel {max(devs):.2e}")


if __name__ == "__main__":
    main()


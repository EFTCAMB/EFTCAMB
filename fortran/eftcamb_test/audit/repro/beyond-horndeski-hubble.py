"""Repro: Beyond-Horndeski model (EFTflag=4, FullMappingEFTmodel=5) expansion history used by CAMB lacks dark energy.

Expected (arXiv:1905.11364 eqs. 3.6, 3.16-3.18): the model's background solver fixes x2 so that
Omega_DE(z=0) = x1+x2+x3+x4 = 1 - Omega_m - Omega_r, hence H(z=0) = H0 (the solver itself reports
'relative error on H_0 ~ 1e-6').  CAMB's H(z) (results.hubble_parameter) should therefore equal H0 at z=0 and stay
close to LCDM (the paper's Fig. 3 shows percent-level changes).
Observed: H(z=0) = H0 sqrt(Omega_m) ~ 0.54 H0, i.e. H(z)^2 = H_LCDM(z)^2 - H0^2 Omega_Lambda: the dark-energy
contribution is missing.  The model's c and Lambda are built as adotoa^2 * f(x_i) from the timestep cache, but the
model does not bind compute_adotoa / compute_dtauda (bindings commented out), so the generic full-model Friedmann
equation is evaluated with adotoa = 0 in the cache.

Run: /Users/marco/Software/miniforge3/bin/python beyond-horndeski-hubble.py
"""

import json
import os
import subprocess
import sys
import tempfile

import numpy as np

TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
        "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
BASE = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/base_params.ini"
COMMON = f"""omnuh2 = 0
massive_neutrinos = 0
massless_neutrinos = 3.046
nu_mass_eigenstates = 0
feedback_level = 0
DEFAULT({BASE})
"""
BH = "EFTflag = 4\nFullMappingEFTmodel = 5\nBeyond_Horndeski_x10 = -0.1\nBeyond_Horndeski_x30 = 1.0\n" \
     "Beyond_Horndeski_x40 = 1.0\nEFT_ghost_stability = T\nEFT_gradient_stability = T\n" + COMMON  # BH2, Table I
GR = "EFTflag = 0\n" + COMMON

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_background(cp)
z = np.array([0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0])
print('JSON' + json.dumps({'z': z.tolist(), 'H': np.asarray(res.hubble_parameter(z)).tolist(), 'H0': cp.H0,
                           'OL': float(res.get_Omega('de', z=0))}))
"""


def run(text: str) -> dict:
    fd, ini = tempfile.mkstemp(suffix=".ini")
    os.close(fd)
    with open(ini, "w") as f:
        f.write(text)
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", WORKER, TREE, ini], capture_output=True, text=True, env=env,
                       timeout=900)
    os.remove(ini)
    js = [ln for ln in p.stdout.splitlines() if ln.startswith("JSON")]
    if not js:
        print(p.stdout[-2000:], p.stderr[-2000:])
        raise SystemExit("CAMB worker failed")
    return {k: np.array(v) if isinstance(v, list) else v for k, v in json.loads(js[0][4:]).items()}


def main() -> None:
    bh, gr = run(BH), run(GR)
    H0, OL = gr["H0"], gr["OL"]
    noDE = np.sqrt(gr["H"] ** 2 - H0**2 * OL)  # LCDM with the Lambda term removed
    print(f"input H0 = {H0}")
    print(f"{'z':>6} {'H_BH(code)':>12} {'H_LCDM':>12} {'H_LCDM w/o DE':>14}")
    for i, z in enumerate(bh["z"]):
        print(f"{z:6.1f} {bh['H'][i]:12.5f} {gr['H'][i]:12.5f} {noDE[i]:14.5f}")
    print(f"observed H_BH(z=0)/H0 = {bh['H'][0] / H0:.5f}   expected 1 (model solver matches H0 to ~1e-6)")
    print(f"max |H_BH / H_noDE - 1| over z<=100: {np.max(np.abs(bh['H'] / noDE - 1)):.2e}  (DE missing if ~0)")
    print("RESULT:", "BUG REPRODUCED" if abs(bh["H"][0] / H0 - 1) > 0.01 else "not reproduced")


if __name__ == "__main__":
    main()

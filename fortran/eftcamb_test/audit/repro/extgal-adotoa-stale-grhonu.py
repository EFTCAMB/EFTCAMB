"""Repro (t2d-static): Extended Galileon expansion rate ignores massive neutrinos for a >= 0.2 (reads a cache field
before any caller has filled it), and its bisection has a zero tolerance (integer 10**(-1) == 0).

Code: fortran/eftcamb/08f_full_models/008p7_Extended_Galileon.f90, EFTCAMBExtendedGalileonComputeAdotoa
  * a <  0.2: adotoa from the generic Friedmann equation with eft_cache%grhom_t (includes massive neutrinos);
  * a >= 0.2: bisection of  Omega_phi0 + Omega_tot (x/a)^S - (x/a)^(2+S) = 0  with
        Omega_tot = (omegac+omegab) a^-3 + (omegag+omegar) a^-4 + eft_cache%grhonu_tot/(3 H0^2 a^2)       (l.455)
    but every caller computes grhonu_tot only AFTER compute_adotoa (equations.f90:2816-2844 derivs, :40-47 dtauda via
    compute_background_EFT_functions -> compute_adotoa (l.338), 09_EFTCAMB_IC.f90:115-140, 11_EFTCAMB_python.f90:
    422-439).  On a freshly initialised cache (dtauda, initial conditions, stability check, get_eft_functions, and every
    derivs call when DebugEFTCAMB) grhonu_tot = 0, i.e. the massive neutrinos are dropped; in release derivs the value
    is whatever the previous derivs call of that k mode left (history dependent).
  * loop condition (l.479)  `sqrt(fsoluction**2) > 10**(-1) .and. counter < 50**1`: 10**(-1) is integer 0, so the loop
    always does 50 halvings of [0, 1e9]: accuracy limited to 1e9/2^51 = 4.4e-7 in a*E.
  * same root cause for a < 0.2 in CAMB's background: dtauda (equations.f90:38-47) passes a fresh cache with only
    grhoa2 set; EFTCAMBExtendedGalileonBackgroundEFTFunctions (l.334-338) calls compute_adotoa, whose a < 0.2 branch
    uses grhom_t, EFTc, EFTLambda (all 0 there) -> adotoa = 0 -> the routine returns before setting any EFT function,
    so dtauda = sqrt(3/grhoa2) has no dark energy for a < 0.2 (column 'dtauda/exact-1', ~ -Omega_DE(a)/2).

Model: fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini (GCCG1, s = 2, q = 0.35) on base_params.ini
(omnuh2 = 0.01, 1 massive eigenstate); stability checks switched off so that the background can be inspected.
Expected: adotoa(a) == exact root of the same equation including massive neutrinos (brentq); continuous at a = 0.2;
adotoa(a=1) == H0.  Control: the same model with massless neutrinos only.

Usage: /Users/marco/Software/miniforge3/bin/python extgal-adotoa-stale-grhonu.py [TREE]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

REPO = Path(__file__).resolve().parents[4]
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0"
)
INI = REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"
A_LIST = [0.1, 0.19, 0.2 - 1e-9, 0.2 + 1e-9, 0.3, 0.5, 0.8, 1.0]

CHILD = r"""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
a = np.array(json.loads(sys.argv[3]))
pars = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_background(pars)
names, arr = pars.EFTCAMB.get_eft_functions(res, a)
pc = res.Params.EFTCAMB_parameter_cache
z = 1 / a - 1
Hz = res.hubble_parameter(z)          # km/s/Mpc, from CAMB's dtauda
out = dict(adotoa=arr["adotoa"].tolist(), grhonu_tot=arr["grhonu_tot"].tolist(),
           adotoa_dtauda=(a * Hz / 299792.458).tolist(),
           pc={k: getattr(pc, k) for k in ("omegac", "omegab", "omegag", "omegar", "omegav", "omegan", "h0_Mpc")})
print("JSON" + json.dumps(out))
"""


def exact_root(pc: dict, a: float, grhonu: float, S: float = 2.0) -> float:
    om_tot = (pc["omegac"] + pc["omegab"]) * a**-3 + (pc["omegag"] + pc["omegar"]) * a**-4 + \
        grhonu / (3 * pc["h0_Mpc"]**2 * a**2)
    x = brentq(lambda x: pc["omegav"] + om_tot * (x / a) ** S - (x / a) ** (2 + S), 1e-6, 1e3, xtol=1e-16, rtol=1e-15)
    return x * pc["h0_Mpc"]


def run(tmp: Path, name: str, extra: str) -> dict:
    w = tmp / f"{name}.ini"
    w.write_text(f"DEFAULT({INI})\nfeedback_level = 0\nEFTCAMB_feedback_level = 0\nEFTCAMB_write_background = F\n"
                 f"EFT_ghost_stability = F\nEFT_gradient_stability = F\nEFT_additional_priors = F\n"
                 f"output_root = {tmp / name}\n{extra}")
    r = subprocess.run([sys.executable, "-c", CHILD, str(TREE), str(w), json.dumps(A_LIST)], capture_output=True,
                       text=True, cwd=tmp, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=600)
    js = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
    if not js:
        print(f"child failed rc={r.returncode}\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
        sys.exit(1)
    return json.loads(js[-1][4:])


def main() -> None:
    tmp = Path(tempfile.mkdtemp(prefix="extgal_"))
    print(f"tree: {TREE}")
    for name, extra in [("massive nu (base_params: omnuh2 = 0.01)", ""),
                        ("control: massless nu only", "omnuh2 = 0\nmassive_neutrinos = 0\nmassless_neutrinos = 3.046\n")]:
        d = run(tmp, name.split()[0], extra)
        pc = d["pc"]
        print(f"== {name}: Omega_nu = {pc['omegan']:.5f}, H0 = {pc['h0_Mpc']:.10e} 1/Mpc")
        print("   a            code adotoa      exact(with nu)   code/exact-1   exact(no nu)/exact-1   dtauda/exact-1")
        for i, a in enumerate(A_LIST):
            ex = exact_root(pc, a, d["grhonu_tot"][i])
            ex0 = exact_root(pc, a, 0.0)
            print(f"   {a:<12.9f} {d['adotoa'][i]:.10e} {ex:.10e} {d['adotoa'][i] / ex - 1:+.3e}     "
                  f"{ex0 / ex - 1:+.3e}             {d['adotoa_dtauda'][i] / ex - 1:+.3e}")
        print(f"   adotoa(a=1)/H0 - 1 = {d['adotoa'][-1] / pc['h0_Mpc'] - 1:+.3e}   (expected 0)")
        print(f"   jump adotoa(0.2-1e-9)/adotoa(0.2+1e-9) - 1 = {d['adotoa'][2] / d['adotoa'][3] - 1:+.3e}"
              f"   (exact: {exact_root(pc, A_LIST[2], d['grhonu_tot'][2]) / exact_root(pc, A_LIST[3], d['grhonu_tot'][3]) - 1:+.3e})")
    print("expected: code/exact-1 ~ 0 (<= 1e-8) at every a, no jump at a = 0.2, adotoa(a=1) = H0;"
          " bisection floor with the zero tolerance: 1e9/2^51/(a E) = 4.4e-7 at a = 1")


if __name__ == "__main__":
    main()

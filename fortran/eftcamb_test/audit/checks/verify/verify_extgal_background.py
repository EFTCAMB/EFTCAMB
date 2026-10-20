"""Verify (adversarial, physics lens) finding extgal-adotoa-stale-grhonu: which code paths and how large.

For the Extended Galileon GCCG1 ini (s=2, q=0.35), trees A0 and FIX3, and two neutrino settings
(base_params omnuh2 = 0.01, and the Planck-default mnu = 0.06 eV), measure against the exact tracker root
E^(2+S) = Omega_phi0 + Omega_tot(a) E^S with Omega_tot including the massive-neutrino density:
  * get_eft_functions adotoa                         (fresh cache path)
  * results.hubble_parameter (CAMB dtauda background)  (fresh cache path)
  * get_redshift_evolution 'H' at k = 0.05/Mpc         (perturbation path: release derivs reuses the cache)
  * derived: thetastar, DAstar, rdrag, age; D_A(z=0.5,1,2)
Usage: python verify_extgal_background.py            (parent: runs both trees, prints a table)
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
SP = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45"
          "/scratchpad/audit")
TREES = {"A0": SP / "trees/A0", "FIX3": SP / "work/t2d-static/FIX3"}
OUT = SP / "work/verify/extgal/background"
INI = REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"
A_LIST = [0.15, 0.2 + 1e-9, 0.3, 0.5, 0.8, 1.0]
NU = {"omnuh2=0.01": "", "mnu=0.06eV": "omnuh2 = 0.000645\n"}

CHILD = r"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, sys.argv[4]); import harness
import numpy as np
with harness.slot(1):
    sys.path.insert(0, sys.argv[1])
    import camb
    a = np.array(json.loads(sys.argv[3]))
    pars = camb.read_ini(sys.argv[2], no_validate=True)
    res = camb.get_background(pars)
    names, arr = pars.EFTCAMB.get_eft_functions(res, a)
    pc = res.Params.EFTCAMB_parameter_cache
    z = 1 / a - 1
    out = dict(adotoa=arr["adotoa"].tolist(), grhonu_tot=arr["grhonu_tot"].tolist(),
               adotoa_dtauda=(a * res.hubble_parameter(z) / 299792.458).tolist(),
               pc={k: getattr(pc, k) for k in ("omegac", "omegab", "omegag", "omegar", "omegav", "omegan", "h0_Mpc")},
               derived={k: v for k, v in res.get_derived_params().items()
                        if k in ("thetastar", "DAstar", "rdrag", "age", "zstar", "rstar")},
               DA=res.angular_diameter_distance(np.array([0.5, 1.0, 2.0])).tolist())
    try:
        ev = res.get_redshift_evolution(0.05, z, ["H"])
        out["adotoa_pert"] = (np.asarray(ev)[:, 0]).tolist()   # CAMB 'H' = adotoa (conformal), 1/Mpc
    except Exception as e:
        out["adotoa_pert_error"] = repr(e)
print("JSON" + json.dumps(out))
"""


def exact_root(pc: dict, a: float, grhonu: float, S: float = 2.0) -> float:
    om_tot = (pc["omegac"] + pc["omegab"]) * a**-3 + (pc["omegag"] + pc["omegar"]) * a**-4 + \
        grhonu / (3 * pc["h0_Mpc"]**2 * a**2)
    x = brentq(lambda x: pc["omegav"] + om_tot * (x / a) ** S - (x / a) ** (2 + S), 1e-6, 1e3, xtol=1e-16, rtol=1e-15)
    return x * pc["h0_Mpc"]


def run(tree: Path, name: str, extra: str) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    w = OUT / f"{name}.ini"
    w.write_text(f"DEFAULT({INI})\nfeedback_level = 0\nEFTCAMB_feedback_level = 0\nEFTCAMB_write_background = F\n"
                 f"EFT_ghost_stability = F\nEFT_gradient_stability = F\nEFT_additional_priors = F\n"
                 f"output_root = {OUT / name}\n{extra}")
    r = subprocess.run([sys.executable, "-c", CHILD, str(tree), str(w), json.dumps(A_LIST),
                        str(REPO / "fortran/eftcamb_test/audit/checks")], capture_output=True, text=True, cwd=OUT,
                       env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900)
    js = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
    if not js:
        raise RuntimeError(f"child failed rc={r.returncode}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    return json.loads(js[-1][4:])


def main() -> None:
    allres = {}
    for nu, extra in NU.items():
        for t, tree in TREES.items():
            d = run(tree, f"{t}_{nu.replace('=', '_')}", extra)
            allres[f"{t}|{nu}"] = d
            pc = d["pc"]
            print(f"== {t}  {nu}: Omega_nu = {pc['omegan']:.5f}")
            print("   a          eftfun/exact-1  dtauda/exact-1  pert(H)/exact-1")
            for i, a in enumerate(A_LIST):
                ex = exact_root(pc, a, d["grhonu_tot"][i])
                pert = d.get("adotoa_pert")
                ps = f"{pert[i] / ex - 1:+.3e}" if pert else d.get("adotoa_pert_error", "n/a")[:60]
                print(f"   {a:<10.6f} {d['adotoa'][i] / ex - 1:+.3e}      {d['adotoa_dtauda'][i] / ex - 1:+.3e}"
                      f"      {ps}")
            print(f"   derived: {d['derived']}  DA(0.5,1,2) = {d['DA']}")
    for nu in NU:
        a, b = allres[f"A0|{nu}"], allres[f"FIX3|{nu}"]
        print(f"== A0/FIX3 - 1, {nu}:")
        for k in a["derived"]:
            print(f"   {k:10s} {a['derived'][k] / b['derived'][k] - 1:+.3e}")
        for zz, x, y in zip((0.5, 1, 2), a["DA"], b["DA"]):
            print(f"   DA(z={zz}) {x / y - 1:+.3e}")
    (OUT / "summary.json").write_text(json.dumps(allres, indent=1))


if __name__ == "__main__":
    main()

"""Adversarial verifier, finding extgal-adotoa-stale-grhonu (independent check).

Different method from the repro:
  * massive-neutrino density taken from CAMB's own background (results.get_background_densities 'nu'), not from the
    EFT cache field grhonu_tot that get_eft_functions fills;
  * CAMB background H(z) (results.hubble_parameter, i.e. dtauda) compared with
      (a) the exact tracker solution with neutrinos,
      (b) a prediction of the stated root cause: dtauda = rho_m(with nu) + rho_DE evaluated with the tracker root
          WITHOUT neutrinos for a >= 0.2, and no DE for a < 0.2.
    If (b) matches to ~1e-6 the attribution (fresh cache: grhonu_tot = 0, grhom_t = EFTc = EFTLambda = 0) is right;
  * H(z=0) vs the input H0; comoving distance to z_star vs the exact integral;
  * the get_scale_evolution (derivs, printing call) adotoa vs exact.
Runs one child per configuration (inside a harness CPU slot, OMP_NUM_THREADS=1).

Usage: python verify_extgal_adotoa_nu.py TREE OUTJSON
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REPO = HERE.parents[4]
INI = REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"
C_KMS = 299792.458

CHILD = r"""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cfg = json.loads(sys.argv[3])
pars = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_background(pars)
a = np.array(cfg["a"])
z = 1 / a - 1
dens = res.get_background_densities(a, vars=["nu", "cdm", "baryon", "photon", "neutrino"])
pc = res.Params.EFTCAMB_parameter_cache
out = dict(H0=pars.H0, H_z0=float(res.hubble_parameter(0.0)), Hz=res.hubble_parameter(z).tolist(),
           dens={k: v.tolist() for k, v in dens.items()},
           pc={k: getattr(pc, k) for k in ("omegac", "omegab", "omegag", "omegar", "omegav", "omegan", "h0_Mpc")})
names, arr = pars.EFTCAMB.get_eft_functions(res, a)
out["adotoa_eftfun"] = arr["adotoa"].tolist()
der = res.get_derived_params()
out["zstar"] = der["zstar"]
out["chistar"] = float(res.comoving_radial_distance(der["zstar"]))
ae = np.array(cfg["a_ev"])
try:
    names, ev = pars.EFTCAMB.get_scale_evolution(res, [0.05], ae)
    out["adotoa_evol"] = np.asarray(ev["adotoa"]).reshape(-1).tolist()
except Exception as e:
    out["adotoa_evol"] = [float("nan")] * len(ae)
    out["evol_error"] = repr(e)
print("JSON" + json.dumps(out))
"""


def tracker(pc: dict, a: float, om_nu: float, S: float = 2.0, with_de: bool = True) -> float:
    """E(a) = H/H0 of E^(2+S) = Omega_phi0 + Omega_tot E^S (the code's equation)."""
    om_tot = (pc["omegac"] + pc["omegab"]) * a**-3 + (pc["omegag"] + pc["omegar"]) * a**-4 + om_nu
    if not with_de:
        return float(np.sqrt(om_tot))
    return brentq(lambda E: pc["omegav"] + om_tot * E**S - E ** (2 + S), 1e-6, 1e9, xtol=1e-15, rtol=1e-15)


def run_child(tree: Path, tmp: Path, name: str, extra: str, a: np.ndarray, a_ev: np.ndarray) -> dict:
    w = tmp / f"{name}.ini"
    w.write_text(f"DEFAULT({INI})\nfeedback_level = 0\nEFTCAMB_feedback_level = 0\nEFTCAMB_write_background = F\n"
                 f"EFT_ghost_stability = F\nEFT_gradient_stability = F\nEFT_additional_priors = F\n"
                 f"output_root = {tmp / name}\n{extra}")
    with harness.slot(1):
        r = subprocess.run([sys.executable, "-c", CHILD, str(tree), str(w),
                            json.dumps({"a": a.tolist(), "a_ev": a_ev.tolist()})], capture_output=True, text=True,
                           cwd=tmp, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900,
                           preexec_fn=lambda: os.nice(10))
    js = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
    if not js:
        raise RuntimeError(f"child failed rc={r.returncode}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    return json.loads(js[-1][4:])


def main() -> None:
    tree = Path(sys.argv[1])
    outjson = Path(sys.argv[2])
    a = np.unique(np.concatenate([np.logspace(-3, 0, 121), [0.15, 0.19, 0.2 - 1e-9, 0.2 + 1e-9, 0.25, 0.3, 0.5,
                                                            0.8, 1.0]]))
    a_ev = np.array([0.1, 0.19, 0.21, 0.3, 0.5, 0.8, 1.0])
    tmp = Path(tempfile.mkdtemp(prefix="vextgal_", dir=outjson.parent))
    configs = {
        "omnuh2_0.01": "",
        "omnuh2_0.00064": "omnuh2 = 0.00064\n",
        "massless": "omnuh2 = 0\nmassive_neutrinos = 0\nmassless_neutrinos = 3.046\n",
    }
    report: dict = {"tree": str(tree)}
    print(f"tree: {tree}")
    for name, extra in configs.items():
        d = run_child(tree, tmp, name, extra, a, a_ev)
        pc = d["pc"]
        h0 = pc["h0_Mpc"]
        H0mpc = d["H0"] / C_KMS
        # Omega_nu(a) as it enters Omega_tot: 8 pi G rho_nu / (3 H0^2), from CAMB's own background
        om_nu = np.array(d["dens"]["nu"]) / a**4 / (3 * h0**2)
        # consistency of the parameter cache with CAMB's background densities (cdm + baryon at a = 1)
        cb = (d["dens"]["cdm"][-1] + d["dens"]["baryon"][-1]) / (3 * h0**2)
        E_ex = np.array([tracker(pc, ai, oi) for ai, oi in zip(a, om_nu)])
        E_nonu = np.array([tracker(pc, ai, 0.0) for ai in a])
        om_tot = (pc["omegac"] + pc["omegab"]) * a**-3 + (pc["omegag"] + pc["omegar"]) * a**-4 + om_nu
        # stated root cause -> predicted dtauda expansion rate
        E_pred = np.where(a >= 0.2, np.sqrt(om_tot + pc["omegav"] * E_nonu**-2.0), np.sqrt(om_tot))
        E_camb = np.array(d["Hz"]) / C_KMS / h0
        E_code = np.array(d["adotoa_eftfun"]) / (a * h0)
        sel = a >= 0.01  # below that DE is < 1e-9 of the total
        r_camb_ex = E_camb / E_ex - 1
        r_camb_pred = E_camb / E_pred - 1
        r_code_ex = E_code / E_ex - 1
        r_code_nonu = E_code / np.where(a >= 0.2, E_nonu, E_ex) - 1
        # comoving distance to z_star: CAMB vs exact tracker with nu (Omega_nu(a) interpolated in log a)
        la = np.log(a)

        def E_exact_of(ai: float) -> float:
            return tracker(pc, ai, float(np.exp(np.interp(np.log(ai), la, np.log(np.maximum(om_nu, 1e-300))))) if
                           np.any(om_nu > 0) else 0.0)

        ast = 1 / (1 + d["zstar"])
        chi_ex = quad(lambda x: 1 / (x * x * E_exact_of(x) * h0), ast, 1.0, limit=400, epsabs=0, epsrel=1e-10)[0]
        E_ev = np.array(d["adotoa_evol"]) / (a_ev * h0)
        E_ev_ex = np.array([E_exact_of(x) for x in a_ev])
        rec = {
            "Omega_nu0": pc["omegan"], "cache_vs_camb_cb": cb / (pc["omegac"] + pc["omegab"]) - 1,
            "H0_input_kms": d["H0"], "H_z0_camb_kms": d["H_z0"], "H_z0/H0-1": d["H_z0"] / d["H0"] - 1,
            "h0_Mpc_vs_H0/c-1": h0 / H0mpc - 1,
            "adotoa_eftfun(a=1)/H0-1": float(E_code[-1] - 1),
            "max|code/exact-1| a>=0.2": float(np.max(np.abs(r_code_ex[a >= 0.2]))),
            "max|code/root_without_nu-1| a>=0.2": float(np.max(np.abs(r_code_nonu[a >= 0.2]))),
            "max|code/exact-1| 0.01<=a<0.2": float(np.max(np.abs(r_code_ex[sel & (a < 0.2)]))),
            "max|dtaudaH/exact-1|": float(np.max(np.abs(r_camb_ex[sel]))),
            "max|dtaudaH/predicted-1| (root-cause model)": float(np.max(np.abs(r_camb_pred[sel]))),
            "chistar_camb": d["chistar"], "chistar_exact": chi_ex, "chistar/exact-1": d["chistar"] / chi_ex - 1,
            "evol_adotoa/exact-1": dict(zip([f"{x:g}" for x in a_ev], (E_ev / E_ev_ex - 1).tolist())),
            "table": {f"{ai:.9f}": {"code/exact-1": float(r_code_ex[i]), "dtauda/exact-1": float(r_camb_ex[i]),
                                    "dtauda/pred-1": float(r_camb_pred[i])}
                      for i, ai in enumerate(a) if ai in (0.1, 0.15, 0.19, 0.2 - 1e-9, 0.2 + 1e-9, 0.3, 0.5, 0.8, 1.0)},
        }
        report[name] = rec
        print(f"== {name}")
        for k, v in rec.items():
            if k != "table":
                print(f"   {k:45s} {v}")
        for k, v in rec["table"].items():
            print(f"   a={k}  code/exact-1={v['code/exact-1']:+.3e}  dtauda/exact-1={v['dtauda/exact-1']:+.3e}"
                  f"  dtauda/rootcause-pred-1={v['dtauda/pred-1']:+.3e}")
    outjson.write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

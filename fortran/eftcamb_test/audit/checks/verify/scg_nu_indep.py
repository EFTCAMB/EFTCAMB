"""Independent check of finding scg-adotoa-drops-massive-nu (verifier).

Different code path from the repro:
  * adotoa is taken from get_scale_evolution (runs the real derivs in equations.f90 with EFTCAMB_printing), not
    from get_eft_functions;
  * CAMB's expansion rate is reconstructed from the conformal-time integral tau(a) (numerical d tau/d a), not from
    hubble_parameter;
  * the massive-neutrino density is taken from results.get_background_densities('nu') (CAMB's own MassiveNu code).
Sharp attribution test: if the root cause is exactly "SCG adotoa omits rho_nu_massive", then
    R(a) = 3 [ (a H_camb)^2 - adotoa_derivs^2 ] / (8 pi G rho_numass a^2)  ==  1  for all a.
Also checks H(z=0) against the input H0.

Usage: python scg_nu_indep.py TREE OUTDIR
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

PARAMS = HERE.parents[2] / "parameters"

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_results(cp)
c = camb.constants.c / 1e3
a = np.logspace(-6, 0, 61)
a = a[a < 0.999999]
a = np.append(a, 1.0)
fields, ev = cp.EFTCAMB.get_scale_evolution(res, [0.01], a)
ad_derivs = np.asarray(ev["adotoa"])[0]
grhom_derivs = np.asarray(ev["grhom_t"])[0]
# CAMB expansion from conformal time: d tau / d a = 1/(a^2 H) => a H = 1/(a dtau/da)
eps = 1e-4
def tau_of_a(x):
    return np.asarray(res.conformal_time(1 / x - 1))
dtda = (tau_of_a(a * (1 + eps)) - tau_of_a(np.minimum(a * (1 - eps), a))) / (a * eps * 2)
# at a=1 use one-sided (z<0 may be fine, but be safe)
dtda[-1] = (tau_of_a(np.array([1 - 1e-4]))[0] - tau_of_a(np.array([1 - 2e-4]))[0]) / 1e-4
aH_tau = 1 / (a * dtda)
aH_hub = np.asarray(res.hubble_parameter(1 / a - 1)) / c * a
dens = res.get_background_densities(a, vars=["photon", "neutrino", "nu", "cdm", "baryon"])
rho_nu_a2 = dens["nu"] / a**2  # 8 pi G rho_numass a^2
rho_rest_a2 = (dens["photon"] + dens["neutrino"] + dens["cdm"] + dens["baryon"]) / a**2
print(json.dumps({
    "a": a.tolist(), "ad_derivs": ad_derivs.tolist(), "grhom_derivs": grhom_derivs.tolist(),
    "aH_tau": aH_tau.tolist(), "aH_hub": aH_hub.tolist(),
    "rho_nu_a2": rho_nu_a2.tolist(), "rho_rest_a2": rho_rest_a2.tolist(),
    "H0_in": cp.H0, "H_z0": float(res.hubble_parameter(0.0)), "omnuh2": cp.omnuh2,
    "num_nu_massive": int(cp.num_nu_massive), "omega_nu": float(res.get_Omega("nu")),
}))
"""


def run(tree: str, extra: dict, outdir: Path, name: str) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    ini = outdir / f"{name}.ini"
    lines = [f"DEFAULT({PARAMS}/5_scg.ini)", f"output_root = {outdir}/{name}", "EFTCAMB_write_background = F",
             "feedback_level = 0", "EFTCAMB_feedback_level = 0"] + [f"{k} = {v}" for k, v in extra.items()]
    ini.write_text("\n".join(lines) + "\n")
    w = outdir / "scg_nu_worker.py"
    w.write_text(WORKER)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        out = subprocess.run([sys.executable, str(w), tree, str(ini)], capture_output=True, text=True, env=env,
                             cwd=outdir, preexec_fn=lambda: os.nice(10), timeout=1200)
    line = [x for x in out.stdout.splitlines() if x.startswith("{")]
    if not line:
        raise RuntimeError(out.stdout[-3000:] + out.stderr[-3000:])
    r = json.loads(line[-1])
    (outdir / f"{name}.json").write_text(json.dumps(r))
    return r


def main() -> None:
    import numpy as np

    tree = sys.argv[1] if len(sys.argv) > 1 else str(harness.TREES["A0"])
    outdir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(tempfile.mkdtemp())
    cases = {
        "shipped_omnuh2_0.01": {},
        "planck_omnuh2_0.00064": {"omnuh2": "0.00064"},
        "massless_only": {"omnuh2": "0", "massive_neutrinos": "0", "massless_neutrinos": "3.046",
                          "nu_mass_eigenstates": "0"},
    }
    summary = {}
    for name, extra in cases.items():
        r = run(tree, extra, outdir, name)
        a = np.array(r["a"])
        ad = np.array(r["ad_derivs"])
        aht = np.array(r["aH_tau"])
        ahh = np.array(r["aH_hub"])
        rnu = np.array(r["rho_nu_a2"])
        print(f"\n== {name}: omnuh2={r['omnuh2']} num_nu_massive={r['num_nu_massive']} "
              f"Omega_nu={r['omega_nu']:.5f}")
        print(f"   H(z=0)/H0_input - 1 = {r['H_z0'] / r['H0_in'] - 1:+.4e};  "
              f"adotoa_derivs(a=1)/(H0/c) - 1 = {ad[-1] / (r['H0_in'] / 2.99792458e5) - 1:+.4e}")
        print("        a     adotoa_derivs/aH_tau-1   aH_hub/aH_tau-1    R=3[(aH)^2-ad^2]/rho_nu a^2")
        rows = []
        for i in range(0, len(a), 6):
            R = 3 * (aht[i] ** 2 - ad[i] ** 2) / (3 * rnu[i]) if rnu[i] > 0 else float("nan")
            # note: (aH)^2 = sum(8piG rho a^2)/3, so the deficit in (aH)^2 should equal rho_nu a^2 / 3
            R = (aht[i] ** 2 - ad[i] ** 2) / (rnu[i] / 3) if rnu[i] > 0 else float("nan")
            rows.append((a[i], ad[i] / aht[i] - 1, ahh[i] / aht[i] - 1, R))
            print(f"   {a[i]:8.2e}   {ad[i] / aht[i] - 1:+.4e}             {ahh[i] / aht[i] - 1:+.3e}"
                  f"          {R:.6f}")
        i = len(a) - 1
        R = (aht[i] ** 2 - ad[i] ** 2) / (rnu[i] / 3) if rnu[i] > 0 else float("nan")
        print(f"   {a[i]:8.2e}   {ad[i] / aht[i] - 1:+.4e}             {ahh[i] / aht[i] - 1:+.3e}          {R:.6f}")
        mask = rnu > 0
        Rall = (aht[mask] ** 2 - ad[mask] ** 2) / (rnu[mask] / 3) if mask.any() else np.array([np.nan])
        summary[name] = {
            "max_abs_rel_adotoa_vs_tau": float(np.max(np.abs(ad / aht - 1))),
            "rel_adotoa_vs_tau_a1e-6": float(ad[0] / aht[0] - 1),
            "rel_adotoa_vs_tau_a1": float(ad[-1] / aht[-1] - 1),
            "Hz0_over_H0_minus1": r["H_z0"] / r["H0_in"] - 1,
            "R_min": float(np.nanmin(Rall)), "R_max": float(np.nanmax(Rall)),
            "max_hub_vs_tau": float(np.max(np.abs(ahh / aht - 1))),
        }
    print("\nSUMMARY", json.dumps(summary, indent=1))
    (outdir / "summary.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

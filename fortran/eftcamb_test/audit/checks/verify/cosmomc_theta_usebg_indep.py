"""Adversarial verification of finding cosmomc-theta-exact-rs-use-background.

Independent method (does not reuse the repro's quadrature):
  * rs_used = cosmomc_theta() * D_M(z*) recovers the sound horizon the Fortran code actually used.
  * Own Simpson integrals on a dense ln(a) grid, with H(a) from hubble_parameter() and the baryon loading taken
    from the CAMBdata fields grhob/grhog (exact R) or 3e4 a ombh2 (CAMB cosmomc definition, approximate R).
  * results.sound_horizon(z*) (table in T, dsound_da_exact Romberg in F) as a second code path.
  * Run for all three 5_hdsk_TG_* models, on A0 (code under audit) and B1 (pre-rebase), to check generality and
    whether the behaviour is inherited.
Also checks the theta -> H0 consequence through set_H0_for_theta(cosmomc_approx=True) on A0 for TG_1.
"""

import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import harness  # noqa: E402

ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
OUT = os.path.join(ROOT, "work", "verify", "cosmomc_theta_usebg")
PARAMS = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"
CHECKS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

WORKER = r"""
import sys, json, os
import numpy as np
from scipy.integrate import simpson
tree, ini, ub, do_h0 = sys.argv[1], sys.argv[2], sys.argv[3] == "T", sys.argv[4] == "1"
sys.path.insert(0, sys.argv[5])
sys.path.insert(0, tree)
import camb
import harness
cp = camb.read_ini(ini, no_validate=True)
keys = harness.parse_ini(ini)
fallback = False
if cp.EFTCAMB.EFTflag == 0 and int(keys.get("EFTflag", "0")) != 0:
    cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)
    fallback = True
cp.EFTCAMB.EFTCAMB_use_background = ub
res = camb.get_background(cp)
ombh2, omdmh2 = cp.ombh2, cp.omch2 + cp.omnuh2
zs = 1048 * (1 + 0.00124 * ombh2**-0.738) * (1 + (0.0783 * ombh2**-0.238 / (1 + 39.5 * ombh2**0.763))
                                             * (omdmh2 + ombh2) ** (0.560 / (1 + 21.1 * ombh2**1.81)))
ast = 1 / (1 + zs)
theta = res.cosmomc_theta()
DM = float(res.angular_diameter_distance(zs)) * (1 + zs)
DMc = float(res.comoving_radial_distance(zs))
rs_used = theta * DM
# own dense quadrature in x = ln a, H in 1/Mpc
x = np.linspace(np.log(1e-9), np.log(ast), 20001)
a = np.exp(x)
H = np.asarray(res.hubble_parameter(1 / a - 1)) / 299792.458
dtau_dx = 1.0 / (a * H)
Rex = 3 * res.grhob * a / (4 * res.grhog)
Rap = 3.0e4 * a * ombh2
r0 = 1.0 / (a[0] * H[0]) / np.sqrt(3.0)  # radiation-era piece below a=1e-9
rs_ex = simpson(dtau_dx / np.sqrt(3 * (1 + Rex)), x=x) + r0
rs_ap = simpson(dtau_dx / np.sqrt(3 * (1 + Rap)), x=x) + r0
out = dict(fallback=fallback, ub=ub, EFTflag=int(cp.EFTCAMB.EFTflag), H0=cp.H0, zs=zs, theta=theta, DM=DM, DMc=DMc,
           rs_used=rs_used, rs_sh=float(res.sound_horizon(zs)), rs_ex=rs_ex, rs_ap=rs_ap,
           Rratio=float(3 * res.grhob / (4 * res.grhog) / (3.0e4 * ombh2)), H1=float(res.hubble_parameter(1.0)))
if do_h0:
    target = 0.010854335
    cp2 = camb.read_ini(ini, no_validate=True)
    if fallback:
        cp2.EFTCAMB.initialize_parameters(cp2, keys, print_header=False)
    cp2.EFTCAMB.EFTCAMB_use_background = ub
    cp2.set_H0_for_theta(target, cosmomc_approx=True, theta_H0_range=(65, 78))
    out["H0_for_theta_target"] = cp2.H0
print("@@" + json.dumps(out))
"""


def run(tree, model, ub, do_h0=False):
    ini = os.path.join(OUT, f"{os.path.basename(tree)}_{model}_ub{ub}.ini")
    with open(ini, "w") as f:
        f.write(f"DEFAULT({PARAMS}/{model}.ini)\nEFTCAMB_use_background = {ub}\nfeedback_level = 0\n")
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-c", WORKER, tree, ini, ub, "1" if do_h0 else "0", CHECKS],
                           capture_output=True, text=True, env=env, cwd=OUT, timeout=900,
                           preexec_fn=lambda: os.nice(10))
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@")]
    if not line:
        return {"error": (p.stdout[-1500:] + p.stderr[-1500:])}
    return json.loads(line[0][2:])


def main():
    os.makedirs(OUT, exist_ok=True)
    trees = sys.argv[1:] or ["A0", "B1"]
    allres = {}
    for t in trees:
        tree = os.path.join(ROOT, "trees", t)
        for model in ("5_hdsk_TG_1", "5_hdsk_TG_2", "5_hdsk_TG_3"):
            r = {ub: run(tree, model, ub, do_h0=(t == "A0" and model == "5_hdsk_TG_1")) for ub in ("F", "T")}
            allres[f"{t}/{model}"] = r
            if "error" in r["F"] or "error" in r["T"]:
                print(t, model, "ERROR", r)
                continue
            F, T = r["F"], r["T"]
            print(f"== {t} {model}  EFTflag={T['EFTflag']} fallback={T['fallback']} H0={T['H0']:.5f}")
            print(f"   H(z=1) F/T rel diff       : {T['H1'] / F['H1'] - 1:+.2e}")
            print(f"   D_M(z*) F/T rel diff      : {T['DM'] / F['DM'] - 1:+.2e}   (D_M vs chi, T: {T['DM'] / T['DMc'] - 1:+.1e})")
            for ub, d in (("F", F), ("T", T)):
                print(f"   ub={ub}: theta={d['theta']:.9f} rs_used={d['rs_used']:.5f}  own rs_approxR={d['rs_ap']:.5f} "
                      f"own rs_exactR={d['rs_ex']:.5f}  sound_horizon(z*)={d['rs_sh']:.5f}")
                print(f"          rs_used/rs_approxR-1={d['rs_used'] / d['rs_ap'] - 1:+.2e}  "
                      f"rs_used/rs_exactR-1={d['rs_used'] / d['rs_ex'] - 1:+.2e}  "
                      f"sound_horizon/rs_exactR-1={d['rs_sh'] / d['rs_ex'] - 1:+.2e}")
            print(f"   theta(T)/theta(F)-1       : {T['theta'] / F['theta'] - 1:+.4e}   "
                  f"own rs_exactR/rs_approxR-1 : {F['rs_ex'] / F['rs_ap'] - 1:+.4e}   R_exact/R_approx={F['Rratio']:.5f}")
            if "H0_for_theta_target" in F:
                print(f"   H0 for cosmomc_theta=0.010854335: F={F['H0_for_theta_target']:.4f}  T={T['H0_for_theta_target']:.4f}"
                      f"  diff={T['H0_for_theta_target'] - F['H0_for_theta_target']:+.4f} km/s/Mpc")
    with open(os.path.join(OUT, "results_" + "_".join(trees) + ".json"), "w") as f:
        json.dump(allres, f, indent=1)


if __name__ == "__main__":
    main()

"""Repro: EFT_IC_type = 2/3 set pi, pi' from a radiation-era super-horizon series (k^2 tau^3) at the EFT turn-on time.

09_EFTCAMB_IC.f90 (cases 2 and 3) evaluates y(EV%w_ix) = h0 (1/4)(1+Omega)/(5+3 Omega) k^2 tau^3 (case 2) at the
perturbation turn-on time (default a = 0.01, matter era, k tau >> 1).  The result is orders of magnitude away from
the quasi-static attractor that EFT_IC_type = 1 uses, so the pi field starts with a huge transient.

The script runs, on the frozen A0 build:
  1) 2_PEFT_Omega_const_2 (Omega = 0.1) with EFT_IC_type = 1 and 2, prints pi just after turn-on at k = 0.01, 0.1, 1
     and the analytic case-2 value, then the change of the CMB / P(k) outputs between IC types 1 and 2;
  2) 2_PEFT_gamma3_exp_1 with EFT_IC_type = 3, which aborts with a RungeKuttaDP45 error.

Run: /Users/marco/Software/miniforge3/bin/python ic23-pi-radiation-series-at-turn-on.py
"""

import os
import subprocess
import sys
import tempfile

ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREE = os.environ.get("REPRO_TREE", ROOT + "/trees/A0")
PARAMS = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"

CHILD = r"""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
ini, ic, ks = sys.argv[2], sys.argv[3], [float(x) for x in sys.argv[4].split(',')]
cp = camb.read_ini(ini, no_validate=True)
out = {}
try:
    res = camb.get_results(cp)
except Exception as e:
    print(json.dumps({'error': repr(e)[:300]})); sys.exit(0)
cl = res.get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)['total'][:2501]
k, z, pk = res.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
a = np.logspace(-2.0, -1.5, 40)
fields, arr = cp.EFTCAMB.get_scale_evolution(res, ks, a)
out = {'cl': cl.tolist(), 'pk': pk[0].tolist(), 'sigma8': float(res.get_sigma8()[0]),
       'a': arr['a'][0].tolist(), 'tau': arr['tau'][0].tolist(), 'pi': arr['pi'].tolist(),
       'Omega': arr['EFTOmegaV'][0].tolist(), 'h0': cp.H0 / 299792.458}
print(json.dumps(out))
"""


def run(model, ic, ks="0.01,0.1,1"):
    tmp = tempfile.mkdtemp()
    ini = os.path.join(tmp, "w.ini")
    with open(ini, "w") as f:
        f.write(f"DEFAULT({PARAMS}/{model}.ini)\noutput_root = {tmp}/r\nfeedback_level = 0\n"
                f"EFTCAMB_write_background = F\nEFT_IC_type = {ic}\n")
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", CHILD, TREE, ini, str(ic), ks], capture_output=True, text=True,
                       cwd=tmp, env=env)
    import json
    for line in p.stdout.splitlines()[::-1]:
        if line.startswith("{"):
            return json.loads(line)
    return {"error": f"rc={p.returncode} " + (p.stderr or p.stdout)[-400:]}


def main():
    import numpy as np

    r1, r2 = run("2_PEFT_Omega_const_2", 1), run("2_PEFT_Omega_const_2", 2)
    ks = [0.01, 0.1, 1.0]
    print("2_PEFT_Omega_const_2 (Omega = 0.1, turn-on a = 0.01)")
    for name, r in (("IC1", r1), ("IC2", r2)):
        if "error" in r:
            print(name, "ERROR", r["error"])
    if "error" not in r1 and "error" not in r2:
        a = np.array(r1["a"])
        i = int(np.argmax(a > 0.0102))
        for ik, k in enumerate(ks):
            tau, om = r2["tau"][i], r2["Omega"][i]
            series = r2["h0"] * 0.25 * (1 + om) / (5 + 3 * om) * k**2 * tau**3
            print(f"  k={k:5.2f}  a={a[i]:.4f}  k*tau={k * tau:8.1f}  pi(IC1)={r1['pi'][ik][i]: .3e}  "
                  f"pi(IC2)={r2['pi'][ik][i]: .3e}  case-2 series value at turn-on={series: .3e}")
        cl1, cl2 = np.array(r1["cl"]), np.array(r2["cl"])
        dtt = np.max(np.abs(cl2[2:, 0] / cl1[2:, 0] - 1))
        dee = np.max(np.abs(cl2[2:, 1] / cl1[2:, 1] - 1))
        dpk = np.max(np.abs(np.array(r2["pk"]) / np.array(r1["pk"]) - 1))
        print(f"  IC2 vs IC1: max|dTT/TT|={dtt:.3e} max|dEE/EE|={dee:.3e} max|dP/P|={dpk:.3e} "
              f"sigma8 {r1['sigma8']:.5f} -> {r2['sigma8']:.5f}")
        print("  expected: pi(IC2) close to the attractor pi(IC1) (|ratio| ~ 1) and spectra changes << 1")
    r3 = run("2_PEFT_gamma3_exp_1", 3, "0.1")
    print("2_PEFT_gamma3_exp_1 with EFT_IC_type = 3:", r3.get("error", "ran OK")[:260])
    print("  expected: a normal run (EFT_IC_type = 1 runs fine for this model)")


if __name__ == "__main__":
    main()

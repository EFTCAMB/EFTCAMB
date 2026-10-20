"""Repro: designer f(R) (B0 = 1e-4) matter power at k > 0.3/Mpc depends at the 1-10% level on EFTCAMB_GR_threshold.

The return-to-GR module (09_EFTCAMB_RGR.f90) moves the perturbation turn-on to the first a where any |EFT function|
exceeds EFTCAMB_GR_threshold (default 1e-8).  For f(R) the modification of gravity on scales below the Compton
wavelength is O(1) (mu -> 4/3) even when Omega ~ 1e-8, so the size of the EFT functions is not a measure of how
GR-like the perturbations are.  The test model 4_Des_F_R_LCDM_2 (B0 = 1e-4) is turned on at a = 0.093 by default;
lowering the threshold (or skipping the RGR) moves the turn-on earlier and changes P(k) at k ~ 1-10 /Mpc by up to
~10%, while the CMB TT/EE change < 1e-4.

Runs the frozen A0 build (3 runs, ~1-3 min each) and prints the RGR time, sigma8 and P(k, z=0) ratios.
Run: /Users/marco/Software/miniforge3/bin/python fr-rgr-threshold-high-k-pk.py
"""

import json
import os
import subprocess
import sys
import tempfile

ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREE = os.environ.get("REPRO_TREE", ROOT + "/trees/A0")
PARAMS = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"
CHILD = r"""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_results(cp)
k, z, pk = res.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
cl = res.get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)
print('JSON' + json.dumps({'a_on': float(res.Params.EFTCAMB.EFTCAMB_pert_turn_on), 'k': k.tolist(), 'pk': pk[0].tolist(),
      'sigma8': float(res.get_sigma8()[0]), 'tt': cl['total'][2:2001, 0].tolist(), 'pp': cl['lens_potential'][2:2001, 0].tolist()}))
"""


def run(extra):
    tmp = tempfile.mkdtemp()
    ini = os.path.join(tmp, "w.ini")
    with open(ini, "w") as f:
        f.write(f"DEFAULT({PARAMS}/4_Des_F_R_LCDM_2.ini)\noutput_root = {tmp}/r\nfeedback_level = 0\n"
                "EFTCAMB_write_background = F\n" + "".join(f"{k} = {v}\n" for k, v in extra.items()))
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", CHILD, TREE, ini], capture_output=True, text=True, cwd=tmp, env=env)
    for line in p.stdout.splitlines():
        if line.startswith("JSON"):
            return json.loads(line[4:])
    raise RuntimeError(p.stderr[-800:])


def main():
    import numpy as np

    cases = [("default (threshold 1e-8)", {}), ("threshold 1e-10", {"EFTCAMB_GR_threshold": "1e-10"}),
             ("threshold 1e-12", {"EFTCAMB_GR_threshold": "1e-12"})]
    out = [(name, run(ov)) for name, ov in cases]
    ref = out[0][1]
    k = np.array(ref["k"])
    sel = [int(np.argmin(abs(k - kk))) for kk in (0.1, 0.3, 1.0, 3.0, 10.0)]
    print("4_Des_F_R_LCDM_2 (designer f(R), B0 = 1e-4)")
    for name, r in out:
        rat = np.array(r["pk"]) / np.array(ref["pk"]) - 1
        tt = np.max(np.abs(np.array(r["tt"]) / np.array(ref["tt"]) - 1))
        pp = np.max(np.abs(np.array(r["pp"]) / np.array(ref["pp"]) - 1))
        print(f"  {name:26s} turn-on a={r['a_on']:.4g}  sigma8={r['sigma8']:.6f}  max|dTT/TT|={tt:.1e} "
              f"max|dCpp/Cpp|={pp:.1e}  P/P_default-1: " + " ".join(f"k={k[i]:.2g}:{rat[i]:+.2e}" for i in sel))
    print("  expected: results converged w.r.t. the GR threshold (changes << 1e-3); observed changes at k >= 1/Mpc "
          "are percent level and not monotonic")


if __name__ == "__main__":
    main()

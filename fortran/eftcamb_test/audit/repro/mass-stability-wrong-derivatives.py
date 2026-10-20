"""Repro: the EFT_mass_stability criterion is evaluated with wrong second derivatives of the EFT functions.

EFTStabilityComputation (09_EFTCAMB_stability.f90:572-575) calls, for a > turn-on and EFT_mass_stability = T,
    call input_model%compute_additional_derivs(a, ...)   ! overwrites EFTcdotdot, EFTLambdadotdot, EFTGamma{1,2,3}PP
    call input_model%compute_mass_factors(a, ...)        ! EFT_mu1, EFT_mu2 use EFTcdotdot, EFTGamma{1,2,3}PP
compute_additional_derivs of the designer base class (06p2_abstract_EFTCAMB_designer.f90:61-257) differentiates, with
dfridr, a helper that calls compute_adotoa(a, ..., eft_cache) with the HOST a and the HOST cache (not temp_a/temp_cache)
and builds cP from the host cache, so the differentiated quantities do not depend on the probe point: the derivatives
come out as 0 (EFTcdotdot = -2*adotoa*EFTcdot). Same pattern in the full-mapping base (06p1_abstract_EFTCAMB_full_map.f90).
For Pure EFT / RPH the model itself already computed these derivatives analytically (EFTCAMBPureEFTstdBackgroundEFTFunctions,
EFTCAMBRPHSecondOrderEFTFunctions); compute_additional_derivs overwrites them with the wrong numbers.

The script prints, for two shipped test models that have EFT_mass_stability = T:
  * the model's own analytic values (camb get_eft_functions, which does not call compute_additional_derivs),
  * the values seen by compute_mass_factors (cache returned by the public EFTTestStability) and mu_1/H^2,
  * the stability verdict of the unmodified ini.
Optional 2nd tree (e.g. the FIX1 copy with temp_a/temp_cache used in the helpers) for before/after.

Usage: /Users/marco/Software/miniforge3/bin/python mass-stability-wrong-derivatives.py [TREE [TREE2]]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
SCR = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit")
TREES = [Path(t) for t in sys.argv[1:]] or [SCR / "trees/A0"]
MODELS = ["3_RPH_MK_power_1", "2_PEFT_Omega_const_1"]
A_LIST = [0.0101, 0.1, 0.5]
FIELDS = ["EFTcdotdot", "EFTGamma1PP", "EFTGamma2PP"]

CHILD = r"""
import sys, json, ctypes, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
from camb import eftcamb
from camb.baseconfig import camblib, CAMBError
alist = json.loads(sys.argv[3]); fields = json.loads(sys.argv[4])
out = {}
pars = camb.read_ini(sys.argv[2], no_validate=True)
try:
    camb.get_background(pars); out["verdict"] = "stable"
except CAMBError as e:
    out["verdict"] = str(e).strip().splitlines()[-1]
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.EFTCAMB.EFT_mass_stability = False
res = camb.get_background(pars)
names, arr = pars.EFTCAMB.get_eft_functions(res, np.array(alist))
out["analytic"] = {f: [float(x) for x in arr[f]] for f in fields}
P = res.Params
P.EFTCAMB.EFT_mass_stability = True
P.EFTCAMB.EFTCAMB_feedback_level = 0
f = camblib.__eftcamb_stability_MOD_eftteststability
f.restype = ctypes.c_int
used = {k: [] for k in fields + ["mu1_over_H2"]}
for a in alist:
    c = eftcamb.EFTCAMB_timestep_cache()
    kmax = ctypes.c_double(10.0)
    f(ctypes.byref(ctypes.c_double(a)), ctypes.byref(kmax), ctypes.byref(P.EFTCAMB),
      ctypes.byref(P.EFTCAMB_parameter_cache), ctypes.byref(c))
    for k in fields:
        used[k].append(getattr(c, k))
    used["mu1_over_H2"].append(c.EFT_mu1 * (a / c.adotoa) ** 2)
out["used"] = used
print("JSON" + json.dumps(out))
"""


def main():
    tmp = Path(tempfile.mkdtemp(prefix="massder_"))
    for tree in TREES:
        print(f"tree: {tree}")
        for m in MODELS:
            w = tmp / f"{m}.ini"
            w.write_text(f"DEFAULT({REPO / 'fortran/eftcamb_test/parameters' / (m + '.ini')})\nfeedback_level = 0\n"
                         f"EFTCAMB_feedback_level = 0\nEFTCAMB_write_background = F\noutput_root = {tmp / m}\n")
            r = subprocess.run([sys.executable, "-c", CHILD, str(tree), str(w), json.dumps(A_LIST), json.dumps(FIELDS)],
                               capture_output=True, text=True, cwd=tmp, env=dict(os.environ, OMP_NUM_THREADS="1"),
                               timeout=900)
            js = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
            if not js:
                print(f"  {m}: child failed rc={r.returncode}\n{r.stdout[-800:]}\n{r.stderr[-800:]}")
                continue
            d = json.loads(js[-1][4:])
            print(f"  {m} (ini has EFT_mass_stability = T): verdict of the unmodified ini = {d['verdict']}")
            for i, a in enumerate(A_LIST):
                txt = "   ".join(f"{k}: model {d['analytic'][k][i]:+.4e} / used {d['used'][k][i]:+.4e}" for k in FIELDS)
                print(f"    a = {a:<7} {txt}   mu1/H^2 used = {d['used']['mu1_over_H2'][i]:+.4e}")
    print("expected: 'used' == 'model' (the analytic derivative of the same EFT functions) for every field; "
          "mass check threshold mu/H^2 < -EFT_mass_stability_rate^2 = -100")


if __name__ == "__main__":
    main()

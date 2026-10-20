"""Compare 21cm C_l from A0 and from the diagnostic tree A0_zfix (GR 00-equation z replaced by the true EFT z in the
21cm line_distortions source, see cm21_build_zfix.py), via the pure Python API.

Control: EFTflag=0 must be bit-identical between the two trees (patch inactive).
Measure: EFT model (pure EFT, Omega = 1 const, as 2_PEFT_Omega_const_1) 21cm C_l change = error of the GR substitution.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

SCR = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit")
OUT = SCR / "work" / "verify" / "21cm"
TREES = {"A0": SCR / "trees" / "A0", "A0_zfix": OUT / "A0_zfix"}

CHILD = r"""
import sys, json, numpy as np
tree, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, tree)
import camb
from camb import constants
from camb.sources import GaussianSourceWindow
assert camb.__file__.startswith(tree), camb.__file__
COSMO = dict(H0=70.0, ombh2=0.0226, omch2=0.112, omk=0.0, tau=0.09, As=2.1e-9, ns=0.96, lmax=600)
EFT = dict(EFTflag=1, PureEFTmodel=1, PureEFTHorndeski=False, EFTwDE=0, PureEFTmodelOmega=1,
           PureEFTmodelGamma1=0, PureEFTmodelGamma2=0, PureEFTmodelGamma3=0, PureEFTmodelGamma4=0,
           PureEFTmodelGamma5=0, PureEFTmodelGamma6=0, EFTOmega0=1.0, EFT_ghost_math_stability=False,
           EFT_mass_math_stability=False, EFT_ghost_stability=True, EFT_gradient_stability=True,
           EFT_mass_stability=True, EFT_mass_stability_rate=10.0, EFT_additional_priors=True,
           EFTCAMB_turn_on_time=0.01, EFTCAMB_stability_time=1e-10, EFTCAMB_stability_threshold=0.0)
res = {}
for eftname, eft in (("GR", dict(EFTflag=0)), ("EFT", EFT)):
    for z in (1.0, 30.0):
        for terms in ("basic+dist", "dist_only", "basic_only"):
            cp = camb.set_params(**COSMO, **eft)
            cp.WantTensors = False
            cp.Want_CMB_lensing = False
            cp.SourceWindows = [GaussianSourceWindow(redshift=z, source_type="21cm", sigma=5 / (constants.f_21cm / 1e6))]
            cp.SourceTerms.line_basic = terms != "dist_only"
            cp.SourceTerms.line_distortions = terms != "basic_only"
            r = camb.get_results(cp)
            res[f"{eftname}|z={z:g}|{terms}"] = r.get_source_cls_dict(raw_cl=True)["W1xW1"][:601]
np.savez(out, **res)
"""

if __name__ == "__main__":
    files = {}
    for name, tree in TREES.items():
        out = OUT / f"zfix_cls_{name}.npz"
        with harness.slot(1):
            p = subprocess.run([sys.executable, "-c", CHILD, str(tree), str(out)], capture_output=True, text=True,
                               env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900, cwd=OUT)
        if p.returncode:
            print(name, "failed", p.stderr[-2000:])
            sys.exit(1)
        files[name] = np.load(out)
    a, b = files["A0"], files["A0_zfix"]
    summary = {}
    ls = np.array([10, 30, 100, 300, 500])
    for key in a.files:
        x, y = a[key], b[key]
        summary[key] = {"bit_identical": bool(np.array_equal(x, y)),
                        "max_rel_l2_600": float(np.max(abs(y[2:] / x[2:] - 1))),
                        "ratio_fix_over_A0": {int(l): float(y[l] / x[l]) for l in ls}}
        print(f"{key:28s} identical={summary[key]['bit_identical']!s:5s} max|fix/A0-1|={summary[key]['max_rel_l2_600']:.3e}"
              "  fix/A0 at l=" + " ".join(f"{l}:{y[l] / x[l]:.4f}" for l in ls))
    for z in ("1", "30"):
        g, e0, e1 = a[f"GR|z={z}|basic+dist"], a[f"EFT|z={z}|basic+dist"], b[f"EFT|z={z}|basic+dist"]
        print(f"z={z}: EFT/GR (A0) l=10,100,500: {e0[10] / g[10]:.4f} {e0[100] / g[100]:.4f} {e0[500] / g[500]:.4f}"
              f" | EFT/GR (zfix) {e1[10] / g[10]:.4f} {e1[100] / g[100]:.4f} {e1[500] / g[500]:.4f}")
    (OUT / "cm21_zfix_compare.json").write_text(json.dumps(summary, indent=1))

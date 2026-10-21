"""Verifier for t3-polterddot-curvature-gr: independent code path (camb.set_cosmology, no ini, lensed CMB, mnu=0.06).

Compares A0 (EFTflag=0), B0 (plain CAMB v2) and the diagnostic tree with AUDIT_POLTERK=3 (both Kf(2) restored), all with
Reion.max_redshift=50 so the separately-known reionization default change is removed.
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

DIAG = harness.AUDIT_ROOT / "work" / "verify" / "polterK" / "attrib_tree"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "polterK" / "pydirect.json"

WORKER = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
omk = float(sys.argv[2])
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.0224, omch2=0.12, omk=omk, mnu=0.06, tau=0.055)
pars.InitPower.set_params(As=2.1e-9, ns=0.965)
pars.set_for_lmax(2000, lens_potential_accuracy=1)
pars.Reion.max_redshift = 50.0
if hasattr(pars, "EFTCAMB"):
    assert pars.EFTCAMB.EFTflag == 0
res = camb.get_results(pars)
tot = res.get_total_cls(2000, CMB_unit="muK", raw_cl=True)
unl = res.get_unlensed_scalar_cls(2000, CMB_unit="muK", raw_cl=True)
print(json.dumps({"tot": tot.tolist(), "unl": unl.tolist(), "camb": camb.__file__}))
"""


def run(tree, omk, mask="0"):
    env = dict(os.environ, OMP_NUM_THREADS="1", AUDIT_POLTERK=mask)
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        out = subprocess.run([sys.executable, "-c", WORKER, str(harness.tree_path(tree)), str(omk)], env=env,
                             capture_output=True, text=True, check=True)
    d = json.loads(out.stdout.strip().splitlines()[-1])
    return {k: np.array(v) for k, v in d.items() if k != "camb"}, d["camb"]


def main():
    summary = {}
    for omk in (0.0, -0.05, 0.03):
        a, fa = run("A0", omk)
        b, fb = run("B0", omk)
        p, fp = run(str(DIAG), omk, "3")
        for lab, x in (("A0", a), ("fix(mask3)", p)):
            for kind in ("unl", "tot"):
                tt = np.abs(x[kind][2:, 0] / b[kind][2:, 0] - 1)
                te = np.abs(x[kind][2:, 3] - b[kind][2:, 3]) / np.sqrt(b[kind][2:, 0] * b[kind][2:, 1])
                ee = np.abs(x[kind][2:, 1] / b[kind][2:, 1] - 1)
                key = f"omk={omk:+.2f} {lab:>10s} vs B0 [{kind}]"
                summary[key] = dict(TT_l2=float(tt[0]), TT_max_l_lt30=float(tt[:28].max()),
                                    TT_max_l_ge100=float(tt[98:].max()), TE_max=float(te.max()),
                                    EE_max=float(ee.max()))
                s = summary[key]
                print(f"{key}: TT(l=2) {s['TT_l2']:.3e}  max TT l<30 {s['TT_max_l_lt30']:.3e}  "
                      f"max TT l>=100 {s['TT_max_l_ge100']:.3e}  max TE/sqrt {s['TE_max']:.3e}  EE {s['EE_max']:.3e}")
    print("camb modules:", fa, fb, fp, sep="\n  ")
    OUT.write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

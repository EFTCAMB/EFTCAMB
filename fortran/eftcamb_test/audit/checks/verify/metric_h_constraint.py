"""Independent check for finding evolve-metric-h-inconsistent.

Different method from the repro (which compares spectra): evolve single k modes with get_scale_evolution
(EFTCAMB_quantities_evolution -> derivs) and measure, along the trajectory of the EFTCAMB_evolve_metric_h=T run, how far
the integrated h' = 2kz is from the 00-constraint value that the default path uses
    h'_00 = (4 k etak + 2 dgrho/Meff2)/(2 - alphaB_BS)/adotoa + s00p dphi' + (s00 + s00k k^2) dphi
(alphaB_BS = -2 * cache alphaB), computed here in Python from the cached quantities of the same derivs call.
For the default (F) run the same formula must reproduce the cached z to round-off (sanity check of the formula).

Usage: python metric_h_constraint.py [tree] [ini] [tag]
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

tree = sys.argv[1] if len(sys.argv) > 1 else "A0"
ini = Path(sys.argv[2]) if len(sys.argv) > 2 else harness.PARAMS_DIR / "5_hdsk_TG_1.ini"
tag = sys.argv[3] if len(sys.argv) > 3 else f"{tree}_{ini.stem}"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "metric_h"
OUT.mkdir(parents=True, exist_ok=True)
KS = [1e-4, 1e-3, 3e-3, 1e-2, 0.1, 1.0, 10.0]

runs = {}
for flag in ("F", "T"):
    npz = OUT / f"{tag}_mh{flag}.npz"
    if not npz.exists():
        st = harness.run_model(tree, ini, npz, overrides={"EFTCAMB_evolve_metric_h": flag}, evol_k=KS, threads=1)
        if st["status"] != "ok":
            raise SystemExit(f"{flag}: {st['status']} {st['log_tail'][-5:]}")
    runs[flag] = harness.load(npz)


def zc(d, ik):
    k = KS[ik]
    g = lambda f: d[f"evol_{f}"][ik]  # noqa: E731
    aB = -2.0 * g("alphaB")
    hdot = (4 * k * g("etak") + 2 * g("dgrho") / g("Meff2")) / (2 - aB) / g("adotoa")
    hdot += g("s00p") * g("delta_phidot") + (g("s00") + g("s00k") * k * k) * g("delta_phi")
    return 0.5 * hdot / k


a = runs["F"]["evol_a"][0]  # cache field "a" (same for all k)
res = {}
sel = [0, 50, 100, 150, 200, 250, 280, 299]
print(f"tree={tree} ini={ini.name}")
print("a samples:", " ".join(f"{a[i]:.1e}" for i in sel))
for ik, k in enumerate(KS):
    zF, zT = runs["F"]["evol_z"][ik], runs["T"]["evol_z"][ik]
    sanity = np.nanmax(np.abs(zc(runs["F"], ik) / zF - 1))
    viol = zT / zc(runs["T"], ik) - 1
    rel = zT / zF - 1
    res[k] = {"formula_sanity_F": float(sanity), "viol_T": viol[sel].tolist(), "zT_over_zF": rel[sel].tolist()}
    print(f"k={k:g}: formula sanity on F run max|zc/z-1|={sanity:.1e}")
    print("   T run: z_evolved/z_00constraint - 1 :", " ".join(f"{x:+.2e}" for x in viol[sel]))
    print("   z(T)/z(F) - 1                       :", " ".join(f"{x:+.2e}" for x in rel[sel]))
(OUT / f"{tag}_constraint.json").write_text(json.dumps({"a": a[sel].tolist(), "res": res}, indent=1))

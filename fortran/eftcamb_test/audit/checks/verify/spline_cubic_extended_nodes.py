"""Verify spline-cubic-endpoints-derivatives: separate end-node null value from derivative inconsistency.
Omega = a as flag-9 spline with nodes extending beyond [0,1] (so a=1 is an interior point) vs analytic Omega=a."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness
import numpy as np

OUT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit/work/verify/spline-cubic-endpoints-derivatives")
OUT.mkdir(parents=True, exist_ok=True)
ini = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"
ref = OUT / "ref.npz"
if not ref.exists():
    print(harness.run_model("A0", ini, ref, threads=1)["status"])
for tag, xs in [("ext8", list(np.linspace(-0.4, 1.4, 8))), ("ext16", list(np.linspace(-0.4, 1.4, 16)))]:
    ov = {"PureEFTmodelOmega": "9", "EFTOmega_Spline_Pixels": str(len(xs))}
    ov.update({f"EFTOmegax{i + 1}": repr(float(x)) for i, x in enumerate(xs)})
    ov.update({f"EFTOmegav{i + 1}": repr(float(x)) for i, x in enumerate(xs)})
    f = OUT / f"{tag}.npz"
    s = harness.run_model("A0", ini, f, overrides=ov, threads=1)
    print(tag, s["status"])
    if s["status"] == "ok":
        c = harness.compare(f, ref)
        a, b = harness.load(f), harness.load(ref)
        ag = a["eft_a_grid"]
        for fld in ("EFTOmegaV", "EFTOmegaP", "EFTOmegaPP"):
            k = "eft_" + fld
            if k in a:
                for at in (0.3, 0.9, 0.99, 1.0):
                    j = np.argmin(abs(ag - at))
                    print(f"  {fld} a={ag[j]:.3f} spline {a[k][j]: .5f} analytic {b[k][j]: .5f}")
        print(f"  worst Cl {c['worst_cl']:.4f} worst Pk {c['worst_pk']:.4f} sigma8 {c['sigma8']:.4f}")

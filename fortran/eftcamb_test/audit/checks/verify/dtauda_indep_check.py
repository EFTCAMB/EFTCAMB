"""Independent check of t3-stability-dtauda-without-matter from the python-bridge EFT functions (A0, skip_stability).

Recomputes the argument of EFTCAMBFullModelComputeDtauda with grhoa2=0 (what the stability check sees) and with
grhoa2 = grhom_t*a^2 (what dtauda in equations.f90 uses); reports sign; also checks the bridge grhoa2 power, and
whether the R5 fix changes accepted spectra (A0 vs fixed at OLOmega0=1e-2, 1e-4)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

D = Path(sys.argv[1])
for name in ("Om0", "Om1e-6", "Om1e-4"):
    d = harness.load(D / f"A0_skipstab_{name}.npz")
    a = d["eft_a_grid"]
    a2 = a * a
    den = 1 + d["eft_EFTOmegaV"] + a * d["eft_EFTOmegaP"]
    de = (2 * d["eft_EFTc"] * a2 - d["eft_EFTLambda"] * a2) / 3 / den
    full = (d["eft_grhom_t"] * a2 + 2 * d["eft_EFTc"] * a2 - d["eft_EFTLambda"] * a2) / 3 / den
    neg = a[de < 0]
    print(f"{name}: DE-only arg <0 at {neg.size}/{a.size} grid points, a range [{neg.min() if neg.size else 0:.2e},"
          f"{neg.max() if neg.size else 0:.2e}]; with matter min arg = {full.min():.3e} (negative pts {np.sum(full < 0)})")
    i = np.argmin(abs(a - 1e-6))
    print(f"   a={a[i]:.2e}: bridge grhoa2={d['eft_grhoa2'][i]:.4e}, grhom_t*a^2={d['eft_grhom_t'][i]*a2[i]:.4e}, "
          f"grhom_t/a^4={d['eft_grhom_t'][i]/a[i]**4:.4e}; c={d['eft_EFTc'][i]:.3e}")
    # a^2 H^2 from adotoa vs full arg (should match: dtauda=1/(a^2H) -> 1/dtauda^2 = a^2 adotoa^2 ... in conformal: adotoa^2*a^2)
    j = np.argmin(abs(a - 1e-3))
    print(f"   a={a[j]:.2e}: (a*adotoa)^2={(a[j]*d['eft_adotoa'][j])**2:.6e} vs full arg={full[j]:.6e}, DE-only={de[j]:.3e}")
for name in ("Om1e-2", "Om1e-4"):
    c = harness.compare(harness.load(D / f"A0_{name}.npz"), harness.load(D / f"fixed_{name}.npz"))
    print(name, "A0 vs R5 bit_identical:", c["bit_identical"], "worst_cl", c["worst_cl"])

"""Repro (t3-gr): EFTflag=0 (GR) low-ell TT differs from plain CAMB v2 in curved models.

The EFTCAMB rewrite of polterddot in equations.f90 (output section of derivs, "polterddot from the photon equations")
    pigddot = 2/5 k qgdot - 3/5 k octgdot - dopacity (pig - polter) - opacity (pigdot - polterdot) + 8/15 k sigmadot
is active for every model, including EFTflag=0. CAMB's own photon quadrupole equation is
    pigdot = denlk(2) qg - denlk2(2) octg - ...,  denlk2(2) = 3/5 k Kf(2),  Kf(2) = 1 - 8 K/k^2,
so the rewrite drops the curvature factor Kf(2) on the octupole term (upstream CAMB's formula keeps it).
In flat space the two forms are algebraically identical (A0 and B0 agree to ~1e-14 once the separate reionization
max_redshift change is removed); with omk/=0 the low-ell TT changes.

Both runs set Reion.max_redshift=50 so that only the polterddot change remains (verified: a build of the code under
audit with only the original polterddot restored and FMA contraction disabled is bit-identical to CAMB v2, also at
omk=+-0.05).
Expected: max |dTT/TT| at ell<30 ~1e-13 as in flat space. Observed: printed.
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "t3-gr" / "repro_polter"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"


def job(a):
    tree, omk = a
    ov = {"omk": omk, "num_redshiftwindows": "0", "get_tensor_cls": "F", "l_max_scalar": "1500", "do_lensing": "F"}
    return harness.run_model(tree, GR, OUT / f"{tree}_omk{omk}.npz", overrides=ov,
                             python_set={"Reion.max_redshift": 50.0}, threads=1)


def main():
    omks = ("0", "-0.05", "0.05", "-0.01")
    with ThreadPoolExecutor(4) as ex:
        for st in ex.map(job, [(t, o) for o in omks for t in ("A0", "B0")]):
            if st["status"] != "ok":
                print(st["out"], st["status"])
    print(f"{'omk':>6s} {'max|dTT/TT| l<30':>18s} {'at l':>5s} {'max l>=100':>11s} {'max|dTE|/sqrt(TT EE)':>21s}")
    for o in omks:
        a = harness.load(OUT / f"A0_omk{o}.npz")["cl_unlensed_scalar"]
        b = harness.load(OUT / f"B0_omk{o}.npz")["cl_unlensed_scalar"]
        r = np.abs(a[2:1500, 0] / b[2:1500, 0] - 1)
        te = np.abs(a[2:1500, 3] - b[2:1500, 3]) / np.sqrt(b[2:1500, 0] * b[2:1500, 1])
        print(f"{o:>6s} {r[:28].max():18.3e} {2 + int(np.argmax(r[:28])):5d} {r[98:].max():11.3e} {te.max():21.3e}")
    print("expected: all ~1e-13 (rounding), as for omk=0")


if __name__ == "__main__":
    main()

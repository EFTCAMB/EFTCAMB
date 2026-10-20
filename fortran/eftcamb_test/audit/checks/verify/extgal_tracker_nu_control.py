"""Control (omnuh2 via env OMNUH2) of extgal_tracker_independent.py.
Verifier (independent route) for t3-extgal-gr-point-background.

Runs EFTCAMB Extended Galileon directly from the A0 tree (inside harness.slot) and compares
  (a) CAMB's background H(z) (results.hubble_parameter, i.e. dtauda) and
  (b) get_eft_functions adotoa
with an independent numpy solution of the tracker equation  E^(2+S) = Omega_phi0 + Omega_tot(a) E^S,
where Omega_tot(a) = rho_noDE(a)/rho_crit0 is taken from CAMB's own grho_no_de (includes massive nu).
"""

import os
import sys
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

sys.path.insert(0, str(harness.TREES["A0"]))
import camb  # noqa: E402
import numpy as np  # noqa: E402

INI = REPO / "fortran" / "eftcamb_test" / "audit" / "checks" / "verify" / "extgal_verify.ini"
C = 299792.458


def run(S, q, omnuh2=None):
    extra = "" if omnuh2 is None else f"omnuh2 = {omnuh2}\n"
    tmp = Path(harness.AUDIT_ROOT) / "work" / "verify" / f"extgal_S{S}_q{q}_nu{omnuh2}.ini"
    tmp.write_text(f"ExtendedGalileon_S = {S}\nExtendedGalileon_q = {q}\nEFTCAMB_skip_stability = T\n"
                   f"{extra}DEFAULT({INI})\n")
    cp = camb.read_ini(str(tmp), no_validate=True)
    print("read ok", file=sys.stderr, flush=True)
    res = camb.get_background(cp)
    print("background ok", file=sys.stderr, flush=True)
    return cp, res


def tracker_E(res, S, a):
    h0 = res.Params.H0 / C  # 1/Mpc
    grhocrit = 3 * h0**2
    names = ["cdm", "baryon", "photon", "neutrino", "nu"]
    aa = np.append(a, 1.0)
    d = res.get_background_densities(aa, vars=names)
    rho = sum(d[n] for n in names) / aa**4 / grhocrit  # Omega_tot(a), no DE, incl. massive nu
    rho_tot, om_phi0 = rho[:-1], 1.0 - rho[-1]
    print("Omega_nu_massive(a=1) =", d["nu"][-1] / grhocrit, flush=True)
    out = []
    for i, ai in enumerate(a):
        f = lambda E: om_phi0 + rho_tot[i] * E**S - E ** (2 + S)  # noqa: E731
        lo, hi = 1e-3, 1e8  # f(lo) > 0 > f(hi); plain float bisection to 1e-15 rel (scipy segfaults beside camblib)
        for _ in range(400):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
            if hi - lo < 1e-15 * hi:
                break
        out.append(0.5 * (lo + hi))
    return np.array(out), rho_tot, om_phi0


def main():
    a = np.array([0.05, 0.1, 0.15, 0.19, 0.21, 0.3, 0.5, 0.8, 1.0])
    z = 1 / a - 1
    with harness.slot(1):
        for S, q in ((2.0, 0.35), (0.0, 100.0)):
            cp, res = run(S, q, omnuh2=float(os.environ.get("OMNUH2", "0.01")))
            E, rho_tot, om_phi0 = tracker_E(res, S, a)
            H0 = res.Params.H0
            print("tracker ok", file=sys.stderr, flush=True)
            Hc = np.array([res.hubble_parameter(zi) for zi in z]) / H0
            try:
                _, arr = cp.EFTCAMB.get_eft_functions(res, a)
                ad = np.array([arr[i]["adotoa"] for i in range(len(a))]) / (H0 / C)
            except Exception as e:  # noqa: BLE001
                print("get_eft_functions failed:", e)
                ad = np.full_like(a, np.nan)
            Enode = np.sqrt(rho_tot)
            print(flush=True)
            print(f"S={S} q={q}  Omega_phi0={om_phi0:.5f}  Omega_nu(a=1) via grho_no_de included")
            print("     a      z    tracker E   H_camb/E-1   adotoa/(aE)-1  noDE E/E-1")
            for i in range(len(a)):
                print(f"  {a[i]:.2f} {z[i]:6.2f}  {E[i]:10.5f}  {Hc[i] / E[i] - 1:+.3e}   "
                      f"{ad[i] / (a[i] * E[i]) - 1:+.3e}    {Enode[i] / E[i] - 1:+.3e}")


if __name__ == "__main__":
    main()

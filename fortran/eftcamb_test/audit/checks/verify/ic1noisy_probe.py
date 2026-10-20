"""Verifier (ic1-pidot-transient-noisy-coefficients): probe the pi-equation coefficients that enter the IC type 1
pi' at the turn-on.  The coefficients C, D1, D2, A1, A2, B1, B2 are background quantities (only E depends on the
perturbations); get_eft_functions does not fill them, so they are read from get_scale_evolution of a run whose
turn-on is 0.5 a_on (so that the whole dfridr stencil a_on (1 +- 0.03) is after the turn-on).
  * smoothness of C, D1 on a dense grid around a_on (dx = 5e-4 in ln a)
  * emulation of the Fortran dfridr (same algorithm/constants, h = 0.03 a) on den(a) = C + k^2 EFTpiDfunction
    -> n_IC = d ln den/d ln a; for Gamma != 0 models E' = 0 in the IC so pi'/(aH pi) = -n_IC.  Done with the shipped
    typo (EFTpiDfunction = D1 + k^2 D1) and with D2.  Compared with +-5% / +-20% secant log-slopes.
  * attractor diagnostics at a_on: M/(A aH^2) and B/(A aH)  (M = den, A = A1 + k^2 A2, B = B1 + k^2 B2)
  * the pi, pi' actually set by the code: get_scale_evolution of the real run (turn-on a_on) at a_on (1 + 1e-6),
    pi'/(aH pi) to compare with -n_IC(typo).
Usage: python ic1noisy_probe.py <tree> <wrapper_half.ini> <wrapper_real.ini> <a_on> <out.json>
"""

import json
import sys

import numpy as np

tree, ini_half, ini_real, a_on, out = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
sys.path.insert(0, tree)
import camb  # noqa: E402

KS = [0.01, 0.1, 0.5, 1.0, 5.0, 10.0]
NSTEN = 60
CON = 1.4
FIELDS = ["EFTpiA1", "EFTpiA2", "EFTpiB1", "EFTpiB2", "EFTpiC", "EFTpiD1", "EFTpiD2", "adotoa", "EFTGamma1V", "EFTc",
          "EFTOmegaV"]

# ---- background coefficients from the half-turn-on run
cp = camb.read_ini(ini_half, no_validate=True)
res = camb.get_results(cp)
xg = np.log(a_on) + np.linspace(-0.1, 0.1, 401)
sten = [a_on]
for i in range(NSTEN):
    h = 0.03 * a_on / CON ** i
    sten += [a_on + h, a_on - h]
a_all = np.unique(np.concatenate([np.exp(xg), np.array(sten)]))
_, arr = cp.EFTCAMB.get_scale_evolution(res, [0.1], a_all)
tab = {f: np.asarray(arr[f][0], dtype=float) for f in FIELDS}
lut = {float(a): i for i, a in enumerate(a_all)}


def coef(a, f):
    return tab[f][lut[float(a)]]


def den(a, k, typo):
    d = coef(a, "EFTpiD1") + k * k * (coef(a, "EFTpiD1") if typo else coef(a, "EFTpiD2"))
    return coef(a, "EFTpiC") + k * k * d


def dfridr(f, x, h, ntab=NSTEN, con=CON, safe=2.0):
    """Numerical Recipes dfridr, same algorithm/constants as fortran/eftcamb/02_utilities.f90 (ntab there is 100;
    the emulation stops at NSTEN steps and reports it)."""
    con2 = con * con
    A = np.zeros((ntab, ntab))
    hh = h
    A[0, 0] = (f(x + hh) - f(x - hh)) / (2 * hh)
    err, ans, n = 1e30, 0.0, 1
    for i in range(1, ntab):
        hh /= con
        n += 1
        A[0, i] = (f(x + hh) - f(x - hh)) / (2 * hh)
        fac = con2
        for j in range(1, i + 1):
            A[j, i] = (A[j - 1, i] * fac - A[j - 1, i - 1]) / (fac - 1)
            fac *= con2
            errt = max(abs(A[j, i] - A[j - 1, i]), abs(A[j, i] - A[j - 1, i - 1]))
            if errt <= err:
                err, ans = errt, A[j, i]
        if abs(A[i, i] - A[i - 1, i - 1]) >= safe * err:
            break
    return ans, err, n, list(A[0, : n])


result = {"a_on": a_on, "k": {}}
for f in FIELDS:
    result[f] = float(coef(a_on, f))
H = result["adotoa"]
dense_idx = [lut[float(a)] for a in np.exp(xg)]
result["dense"] = {"x": xg.tolist(), **{f: tab[f][dense_idx].tolist() for f in ("EFTpiC", "EFTpiD1", "EFTpiD2",
                                                                              "EFTpiA1", "adotoa")}}
for k in KS:
    M = result["EFTpiC"] + k * k * result["EFTpiD1"] + k ** 4 * result["EFTpiD2"]
    Aa = result["EFTpiA1"] + k * k * result["EFTpiA2"]
    Bb = result["EFTpiB1"] + k * k * result["EFTpiB2"]
    row = {"M": M, "M_over_A_H2": M / Aa / H ** 2 if Aa else None, "B_over_AH": Bb / Aa / H if Aa else None}
    for typo in (True, False):
        tag = "typo" if typo else "D2"
        f = lambda a, k=k, typo=typo: den(a, k, typo)  # noqa: E731
        d, err, n, cdiffs = dfridr(f, a_on, 0.03 * a_on)
        den0 = f(a_on)
        row[f"n_IC_{tag}"] = d * a_on / den0
        row[f"dfridr_relerr_{tag}"] = err / abs(d) if d else None
        row[f"dfridr_nsteps_{tag}"] = n
        row[f"centered_logslopes_{tag}"] = [c * a_on / den0 for c in cdiffs]
        for frac in (0.05, 0.2):
            lo, hi = a_on * np.exp(-frac), a_on * np.exp(frac)
            ilo, ihi = np.argmin(abs(np.exp(xg) - lo)), np.argmin(abs(np.exp(xg) - hi))
            dl, dh = (tab["EFTpiC"][dense_idx[i]] + k * k * (tab["EFTpiD1"][dense_idx[i]] + k * k * (
                tab["EFTpiD1"][dense_idx[i]] if typo else tab["EFTpiD2"][dense_idx[i]])) for i in (ilo, ihi))
            row[f"n_secant{int(frac * 100)}_{tag}"] = float((np.log(abs(dh)) - np.log(abs(dl))) / (xg[ihi] - xg[ilo])
                                                            ) if dl * dh > 0 else "sign change"
    result["k"][str(k)] = row

# ---- the IC actually set by the code (real turn-on)
cp2 = camb.read_ini(ini_real, no_validate=True)
res2 = camb.get_results(cp2)
a_ic = np.array([a_on * (1 + 1e-6), a_on * 1.01, a_on * 1.1, a_on * 2.0, a_on * 10.0])
_, arr2 = cp2.EFTCAMB.get_scale_evolution(res2, KS, a_ic)
result["ic_actual"] = {}
for ik, k in enumerate(KS):
    pi, pid, Hh = (np.asarray(arr2[f][ik], dtype=float) for f in ("pi", "pidot", "adotoa"))
    M = np.asarray(arr2["EFTpiC"][ik]) + k * k * np.asarray(arr2["EFTpiD1"][ik]) + k ** 4 * np.asarray(
        arr2["EFTpiD2"][ik])
    piqs = -cp2.H0 / 299792.458 * np.asarray(arr2["EFTpiE"][ik]) / M
    result["ic_actual"][str(k)] = {"a": a_ic.tolist(), "pi": pi.tolist(), "pidot": pid.tolist(),
                                   "pidot_over_aHpi": (pid / Hh / pi).tolist(), "pi_QS": piqs.tolist()}
json.dump(result, open(out, "w"))
print(json.dumps({kk: v for kk, v in result.items() if kk != "dense"}, indent=1)[:6000])

"""Independent check of the per-model derivative catalogue in finding fullmap-higher-derivs-unset.

Method (different from the t4b-derivs scripts and the repro): read the reference harness outputs (runs/A0/*.npz, the
eft_* timestep-cache fields returned by get_eft_functions on a = logspace(-8, 0, 401)) and, for every model, compare
each stored derivative with a cubic-spline derivative (in ln a) of the next-lower stored field:

    Omega, Gamma1     : XP  = (1/a) dX/dln a
    c, Lambda         : Xdot = adotoa dX/dln a - 2 adotoa X   (EFTCAMB convention X = x a^2, Xdot = a^2 x')
                        (convention verified on the control pair EFTcdot <- EFTc and on Pure-EFT models,
                         where EFTcdotdot / EFTLambdadotdot are analytic)
    adotoa, Hdot, ... : d/dtau = adotoa d/dln a

Only first derivatives of splines are used.  Errors are max |stored - spline| / scale over 1e-2 <= a <= 0.95, with
scale = max |spline derivative| (and, for H, max(|FD|, adotoa^{n+1}) pointwise as in the finding).  A field is flagged
"ZERO-BUT-NEEDED" when it is identically 0 while its spline derivative is non-negligible relative to the base field.

Usage: python hderivs_catalogue.py [runs_dir] > out.txt
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

RUNS = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/runs/A0")
AMIN, AMAX = 1e-2, 0.95

PAIRS = [  # (base, deriv, kind)
    ("EFTOmegaV", "EFTOmegaP", "a"),
    ("EFTOmegaP", "EFTOmegaPP", "a"),
    ("EFTOmegaPP", "EFTOmegaPPP", "a"),
    ("EFTOmegaPPP", "EFTOmegaPPPP", "a"),
    ("EFTGamma1V", "EFTGamma1P", "a"),
    ("EFTGamma1P", "EFTGamma1PP", "a"),
    ("EFTc", "EFTcdot", "a2"),
    ("EFTcdot", "EFTcdotdot", "a2"),
    ("EFTLambda", "EFTLambdadot", "a2"),
    ("EFTLambdadot", "EFTLambdadotdot", "a2"),
    ("adotoa", "Hdot", "tau"),
    ("Hdot", "Hdotdot", "tau"),
    ("Hdotdot", "Hdotdotdot", "tau"),
]
HORDER = {"Hdot": 2, "Hdotdot": 3, "Hdotdotdot": 4}


def spline_deriv(x, y, lna):
    return CubicSpline(lna, y)(lna, 1)


def analyse(npz):
    d = np.load(npz, allow_pickle=True)
    meta = json.loads(str(d["meta"]))
    a = d["eft_a_grid"]
    lna = np.log(a)
    H = d["eft_adotoa"]
    sel = (a >= AMIN) & (a <= AMAX)
    out = {"model": meta.get("eft_model_name"), "EFTflag": meta.get("EFTflag"), "pairs": {}}
    for base, der, kind in PAIRS:
        X = d[f"eft_{base}"]
        Y = d[f"eft_{der}"]
        dX = spline_deriv(a, X, lna)
        if kind == "a":
            fd = dX / a
        elif kind == "a2":
            fd = H * dX - 2 * H * X
        else:
            fd = H * dX
        fds, Ys = fd[sel], Y[sel]
        if kind == "tau":
            scale = np.maximum(np.abs(fds), H[sel] ** HORDER[der])
            err = float(np.max(np.abs(Ys - fds) / scale))
        else:
            sc = np.max(np.abs(fds))
            err = float(np.max(np.abs(Ys - fds)) / sc) if sc > 0 else (0.0 if not np.any(Ys) else np.inf)
        # significance of the base-field variation: |FD| relative to the base field magnitude (per unit ln a)
        if kind == "a":
            signif = float(np.max(np.abs(fds * a[sel])) / max(np.max(np.abs(X[sel])), 1e-300))
            base_mag = float(np.max(np.abs(X[sel])))
        elif kind == "a2":
            signif = float(np.max(np.abs(fds)) / max(np.max(np.abs(H[sel] * X[sel])), 1e-300))
            base_mag = float(np.max(np.abs(X[sel])))
        else:
            signif = 1.0
            base_mag = float(np.max(np.abs(X[sel])))
        allzero = not np.any(Y[sel] != 0)
        out["pairs"][der] = {"err": err, "allzero": bool(allzero), "base_mag": base_mag, "signif": signif,
                             "max|stored|": float(np.max(np.abs(Ys))), "max|fd|": float(np.max(np.abs(fds)))}
    return out


def main():
    files = sorted(RUNS.glob("*.npz"))
    fields = [p[1] for p in PAIRS]
    res = {}
    for f in files:
        try:
            res[f.stem] = analyse(f)
        except Exception as e:  # noqa: BLE001
            print(f"# {f.stem}: {e}")
    short = {"EFTOmegaP": "OmP", "EFTOmegaPP": "OmPP", "EFTOmegaPPP": "OmPPP", "EFTOmegaPPPP": "OmP4",
             "EFTGamma1P": "G1P", "EFTGamma1PP": "G1PP", "EFTcdot": "cd", "EFTcdotdot": "cdd",
             "EFTLambdadot": "Ld", "EFTLambdadotdot": "Ldd", "Hdot": "Hd", "Hdotdot": "Hdd", "Hdotdotdot": "Hddd"}
    print("rel. error of stored derivative vs spline derivative of the stored lower field, 1e-2<=a<=0.95")
    print("  '0!' = stored identically 0 although base varies (signif > 1e-3);  '0 ' = stored 0 and base ~constant/0")
    print(f"{'model':32s} {'name':22s} " + " ".join(f"{short[x]:>9s}" for x in fields))
    for m, r in res.items():
        cells = []
        for x in fields:
            p = r["pairs"][x]
            if p["allzero"]:
                if p["base_mag"] == 0 or p["signif"] < 1e-3:
                    cells.append(f"{'0':>9s}")
                else:
                    cells.append(f"{'0!' + format(p['err'], '.0e'):>9s}")
            else:
                cells.append(f"{p['err']:9.1e}")
        print(f"{m:32s} {str(r['model'])[:22]:22s} " + " ".join(cells))
    Path(sys.argv[2] if len(sys.argv) > 2 else "/dev/null").write_text(json.dumps(res, indent=1)) if len(
        sys.argv) > 2 else None


if __name__ == "__main__":
    main()

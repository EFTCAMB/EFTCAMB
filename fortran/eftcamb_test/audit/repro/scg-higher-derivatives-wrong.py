"""Repro (track t4b-derivs): Scaling Cubic Galileon (EFTflag=4, FullMappingEFTmodel=6) third derivatives are wrong.

fortran/eftcamb/08f_full_models/008p6_Scaling_Cubic.f90, EFTCAMBScalingCubicBackgroundEFTFunctions /
SecondOrderEFTFunctions / ComputeHubbleDer:
  * hpprime (= d^2 h/dN^2, N = ln a, Hdot = H^2 (1+h)) has the wrong overall sign and the wrong sign of 2 A lambda in
    the y2 coefficient (1 + 2 A lambda + 2 A beta2 instead of 1 - 2 A lambda + 2 A beta2), in both copies;
  * fpprime misses the (2 A lambda - 1)(sqrt6 - 6 A x1) x1pp term, so x1ppp is wrong;
  * EFTLambdadotdot is not the derivative of EFTLambdadot (stray constant 3(2 A lambda - 1), stray factor a, hpprime
    multiplied by x1p instead of x1, ...);
  * EFTGamma2PPP has 3*h*2 where 3*h**2 is needed (and uses the wrong hpprime/x1ppp).
Hence Hdotdotdot, EFTLambdadotdot and EFTGamma2PPP returned by get_eft_functions (and used by the positivity-bound
prior, compute_positivity_bounds in 06_abstract_EFTCAMB_model.f90) are not the derivatives of Hdotdot,
EFTLambdadot and EFTGamma2PP.  First and second derivatives (Hdot, Hdotdot, EFTLambdadot, EFTGamma2P, EFTGamma2PP)
are shown as controls: they agree with the numerical derivative.

Method: get_eft_functions on a 9-point stencil in ln a around each a (steps h and h/2, 8th-order central
differences, Richardson-combined); conventions: Xdot = d/dtau = adotoa d/dN for H derivatives,
EFTLambdadotdot = adotoa (dEFTLambdadot/dN - 2 EFTLambdadot), EFTGamma2P.. = d/da = (1/a) d/dN.
An independent symbolic check of every hand-written formula is fortran/eftcamb_test/audit/checks/t4b-derivs/
scg_formula_sympy.py.

Usage: /Users/marco/Software/miniforge3/bin/python scg-higher-derivatives-wrong.py [TREE]
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

try:  # share the audit machine-wide CPU slots when the harness is available (falls back to no limit)
    sys.path.insert(1, str(Path(__file__).resolve().parents[1] / "checks"))
    from harness import slot
except Exception:  # noqa: BLE001
    from contextlib import nullcontext as slot

TREE = sys.argv[1] if len(sys.argv) > 1 else (
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
    "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
PARAMS = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "parameters"))

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_background(cp)
ac = np.array([1e-3, 1e-2, 0.05, 0.2, 0.5, 0.8])
offs = np.arange(-4, 5)
hs = np.array([0.01, 0.005])
x = np.log(ac)[None, None, :] + hs[:, None, None] * offs[None, :, None]
a_all = np.exp(x).ravel()
a_u, inv = np.unique(a_all, return_inverse=True)
_, arr = cp.EFTCAMB.get_eft_functions(res, a_u)
out = {"ac": ac.tolist(), "hs": hs.tolist()}
for f in ["adotoa", "Hdot", "Hdotdot", "Hdotdotdot", "EFTLambda", "EFTLambdadot", "EFTLambdadotdot",
          "EFTGamma2V", "EFTGamma2P", "EFTGamma2PP", "EFTGamma2PPP"]:
    out[f] = np.asarray(arr[f]).ravel()[inv].reshape(x.shape).tolist()
print(json.dumps(out))
"""

W8 = [1 / 280, -4 / 105, 1 / 5, -4 / 5, 0.0, 4 / 5, -1 / 5, 4 / 105, -1 / 280]


def dN(st, hs):
    """st[h][j][c] -> Richardson-combined d/dN at the centres (8th order: (2^8 D(h/2) - D(h)) / (2^8 - 1))."""
    d = [[sum(W8[j] * st[ih][j][c] for j in range(9)) / hs[ih] for c in range(len(st[0][0]))] for ih in range(2)]
    return [(256 * d[1][c] - d[0][c]) / 255 for c in range(len(d[0]))], d


def main() -> None:
    tmp = tempfile.mkdtemp(prefix="scg_hddd_")
    ini = os.path.join(tmp, "scg.ini")
    with open(ini, "w") as f:
        f.write(f"DEFAULT({PARAMS}/5_scg.ini)\noutput_root = {tmp}/out\nEFTCAMB_write_background = F\n"
                "feedback_level = 0\nEFTCAMB_feedback_level = 0\n")
    w = os.path.join(tmp, "w.py")
    with open(w, "w") as f:
        f.write(WORKER)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    with slot(1):
        p = subprocess.run([sys.executable, w, TREE, ini], capture_output=True, text=True, env=env, cwd=tmp,
                           preexec_fn=lambda: os.nice(10))
    line = [x for x in p.stdout.splitlines() if x.startswith("{")]
    if not line:
        raise RuntimeError(p.stdout[-2000:] + p.stderr[-2000:])
    r = json.loads(line[-1])
    ac, hs = r["ac"], r["hs"]
    c0 = lambda f: [r[f][0][4][i] for i in range(len(ac))]  # noqa: E731
    H = c0("adotoa")
    rows = [  # (stored, base, convention)
        ("Hdot", "adotoa", "tau"), ("Hdotdot", "Hdot", "tau"), ("Hdotdotdot", "Hdotdot", "tau"),
        ("EFTLambdadot", "EFTLambda", "tau_a2"), ("EFTLambdadotdot", "EFTLambdadot", "tau_a2"),
        ("EFTGamma2P", "EFTGamma2V", "a"), ("EFTGamma2PP", "EFTGamma2P", "a"), ("EFTGamma2PPP", "EFTGamma2PP", "a"),
    ]
    print("tree:", TREE)
    print("model: 5_scg.ini (Scaling Cubic Galileon, default parameters)")
    print("rel.err = |stored - numerical| / max(|numerical|, natural scale H|X| or |X|/a);  FD-err = |D(h) - D(h/2)| / "
          "|numerical| (accuracy of the numerical derivative)")
    print(f"{'field':16s} " + " ".join(f"a={a:<9.3g}" for a in ac))
    for stored, base, conv in rows:
        D, Dh = dN(r[base], hs)
        S = c0(stored)
        X = c0(base)
        num = []
        for i, a in enumerate(ac):
            if conv == "tau":
                num.append(H[i] * D[i])
            elif conv == "tau_a2":
                num.append(H[i] * (D[i] - 2 * X[i]))
            else:
                num.append(D[i] / a)
        nat = [H[i] * abs(X[i]) if conv != "a" else abs(X[i]) / ac[i] for i in range(len(ac))]
        rel = [abs(S[i] - num[i]) / max(abs(num[i]), nat[i]) for i in range(len(ac))]
        fde = [abs(Dh[0][i] - Dh[1][i]) / max(abs(D[i]), 1e-300) for i in range(len(ac))]
        print(f"{stored:16s} " + " ".join(f"{e:<11.2e}" for e in rel) + "   (FD-err max "
              f"{max(fde):.1e})")
    print("expected: rel.err ~< 1e-6 for every row (as for the control rows Hdot, Hdotdot, EFTLambdadot, EFTGamma2P(P))")
    i = 1
    print(f"example a={ac[i]}: Hdotdotdot stored={c0('Hdotdotdot')[i]:+.6e}, "
          f"numerical={H[i] * dN(r['Hdotdot'], hs)[0][i]:+.6e}; EFTLambdadotdot stored={c0('EFTLambdadotdot')[i]:+.6e}, "
          f"numerical={H[i] * (dN(r['EFTLambdadot'], hs)[0][i] - 2 * c0('EFTLambdadot')[i]):+.6e}")


if __name__ == "__main__":
    main()

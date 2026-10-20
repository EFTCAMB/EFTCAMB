"""Summarise an evol.py npz: per k, pi growth and oscillation diagnostics (track t10b-convergence).

python evsummary.py EVOL.npz [--fields pi,EFTpiA1,...] [--a 1e-4,1e-3,...]
"""

from __future__ import annotations

import argparse
import json

import numpy as np


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("npz")
    p.add_argument("--fields", default="pi,pidot,EFTpiA1,EFTpiA2,EFTpiB1,EFTpiB2,EFTpiC,EFTpiD1,EFTpiD2,EFT_kinetic,"
                   "EFT_gradient,EFT_mu1,EFT_mu2,EFTISW,Psi,Phi")
    p.add_argument("--a", default="1e-5,1e-4,1e-3,3e-3,1e-2,3e-2,0.1,0.3,0.5,1")
    a = p.parse_args()
    d = np.load(a.npz)
    ks, ag = d["evol_k"], d["evol_a"]
    marks = [float(x) for x in a.a.split(",")]
    idx = [int(np.argmin(np.abs(np.log(ag / m)))) for m in marks]
    np.set_printoptions(precision=3, linewidth=250)
    print("a      ", np.array(ag[idx]))
    for i, k in enumerate(ks):
        pi = d["evol_pi"][i]
        sgn = np.sign(pi)
        nz = int(np.sum(sgn[1:] * sgn[:-1] < 0))
        tau = d["evol_tau"][i]
        print(f"--- k={k:g}  pi sign changes={nz}")
        for f in a.fields.split(","):
            if f"evol_{f}" in d:
                print(f"  {f:13s}", np.array(d[f'evol_{f}'][i][idx]))
        # local oscillation frequency from the pi equation, omega^2 = (C + k^2 D1 + k^4 D2)/A, A = A1 + k^2 A2
        A = d["evol_EFTpiA1"][i] + k**2 * d["evol_EFTpiA2"][i]
        M = d["evol_EFTpiC"][i] + k**2 * d["evol_EFTpiD1"][i] + k**4 * d["evol_EFTpiD2"][i]
        with np.errstate(all="ignore"):
            w2 = M / A
        print("  omega^2*tau^2", np.array((w2 * tau**2)[idx]))
        print("  omega/(aH)   ", np.array((np.sign(w2) * np.sqrt(np.abs(w2)) / d["evol_adotoa"][i])[idx]))


if __name__ == "__main__":
    main()

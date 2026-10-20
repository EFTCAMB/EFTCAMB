"""Repro: Horava model Hdot / Hdotdot are inconsistent with its own adotoa (08f_full_models/008p1_Horava.f90).

compute_adotoa sets adotoa**2 = F (grhom_t + grhov_t)/3 with F = (2+2 xi-eta)/(3 lambda+2), so the consistent
derivative is Hdot = d adotoa/dtau = -(adotoa**2 + F (gpresm_t + gpiv_t))/2.  compute_H_derivs instead computes the GR
expression -(adotoa**2 + gpresm_t + gpiv_t)/2 with the already-rescaled adotoa and multiplies it by sqrt(F).
The second Friedmann equation in EFT form, (1+Omega)(2 Hdot + adotoa^2) + gpresm_t + Lambda = 0, is then violated.

Usage: python horava-hdot-inconsistent.py [tree]   (default: frozen A0 tree)
"""

import os
import sys
import tempfile
from pathlib import Path

import numpy as np

AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else AUDIT / "trees" / "A0"
PARAMS = Path(__file__).resolve().parents[2] / "parameters"
sys.path.insert(0, str(TREE))
os.environ["OMP_NUM_THREADS"] = "1"
try:  # share the audit machine-wide CPU slots when the harness is available (falls back to no limit)
    sys.path.insert(1, str(Path(__file__).resolve().parents[1] / "checks"))
    from harness import slot
except Exception:  # noqa: BLE001
    from contextlib import nullcontext as slot
import camb  # noqa: E402

W8 = np.array([1 / 280, -4 / 105, 1 / 5, -4 / 5, 0.0, 4 / 5, -1 / 5, 4 / 105, -1 / 280])


def run(model: str, a_c: np.ndarray, h: float = 0.01):
    tmp = Path(tempfile.mkdtemp(prefix="t4b_horava_"))
    ini = tmp / "w.ini"
    ini.write_text(f"DEFAULT({PARAMS / model}.ini)\noutput_root = {tmp}/out\nEFTCAMB_write_background = F\n"
                   "EFTCAMB_feedback_level = 0\nfeedback_level = 0\n")
    cp = camb.read_ini(str(ini), no_validate=True)
    res = camb.get_background(cp)
    offs = np.arange(-4, 5)
    a = (a_c[:, None] * np.exp(h * offs[None, :])).ravel()
    _, arr = cp.EFTCAMB.get_eft_functions(res, a)
    g = {f: arr[f].reshape(len(a_c), 9) for f in arr.dtype.names}
    H = g["adotoa"][:, 4]
    Hdot_fd = H * (g["adotoa"] @ W8) / h            # d adotoa / dtau = adotoa d adotoa / d ln a
    Hdotdot_fd = H * (g["Hdot"] @ W8) / h            # d Hdot / dtau  (of the stored Hdot)
    c = {f: g[f][:, 4] for f in g}
    # F as defined by compute_adotoa: adotoa**2 = F (grhom_t + grhov_t)/3  (constant in a)
    F = float(np.median(3 * H**2 / (c["grhom_t"] + c["grhov_t"])))
    Hdot_consistent = -0.5 * (H**2 + F * (c["gpresm_t"] + c["gpiv_t"]))
    fried2 = (1 + c["EFTOmegaV"]) * (2 * c["Hdot"] + H**2) + c["gpresm_t"] + c["EFTLambda"]
    return F, c, H, Hdot_fd, Hdotdot_fd, Hdot_consistent, fried2


def main():
    a_c = np.array([1e-3, 0.01, 0.1, 0.5, 0.9])
    for model in ("5_Horava_SolSyst_3", "5_Horava_eta_4", "5_Horava_eta_1"):
        F, c, H, Hdot_fd, Hdotdot_fd, Hdot_cons, fried2 = run(model, a_c)
        print(f"== {model}: F = (2+2xi-eta)/(3lambda+2) = {F:.6f}   (tree {TREE.name})")
        print("        a      Hdot(stored)    dH/dtau(FD)   -(H^2+F p)/2   rel.err(stored)  2nd-Friedmann/H^2")
        for i, a in enumerate(a_c):
            rel = abs(c["Hdot"][i] - Hdot_fd[i]) / abs(Hdot_fd[i])
            print(f"   {a:8.3g}  {c['Hdot'][i]: .6e}  {Hdot_fd[i]: .6e}  {Hdot_cons[i]: .6e}   {rel:.3e}        "
                  f"{fried2[i] / H[i]**2: .3e}")
        rel2 = np.abs(c["Hdotdot"] - Hdotdot_fd) / np.abs(Hdotdot_fd)
        print(f"   Hdotdot(stored) vs d(Hdot)/dtau: max rel err {rel2.max():.3e};  expected: FD ~ 1e-10, "
              f"2nd Friedmann ~ 1e-15")


if __name__ == "__main__":
    with slot(1):
        main()

"""Verifier for t3-gamma1-mass-instability-near-gr.

Independent evaluation of EFT_mu2 (transcribed from the Fortran under audit) for a constant-Gamma1 Pure EFT model on
an LCDM background (Omega = 0, c = 0, Gamma2..6 = 0), in
  (1) double precision (as the code does) and
  (2) 60-digit mpmath,
with the background taken from the CAMB cache (Gamma1 = 1e-2 run, a ~ 0.01) and from an exact matter+Lambda
background.  If mpmath gives a Gamma1-independent finite mu_2 while double precision blows up like 1/Gamma1^2, the
rejection is a rounding/cancellation artefact; if mpmath also diverges, the formula itself is singular in the GR limit.
"""

import math
import re
import sys

import mpmath as mp
import numpy as np

SRC = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb/06_abstract_EFTCAMB_model.f90"
NPZ = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
       "scratchpad/audit/work/t3-gr/repro_gamma1_mass/G1_0.01.npz")
mp.mp.dps = 60


def transcribe():
    text = open(SRC).read()
    body = text[text.index("subroutine EFTCAMBModelComputeMassFactors"):
                text.index("end subroutine EFTCAMBModelComputeMassFactors")]
    body = body[body.index("rho_d = "):]
    joined = ""
    for ln in body.splitlines():
        s = ln.split("!")[0].strip()
        if not s:
            continue
        joined = joined[:-1] + " " + s if joined.endswith("&") else joined + "\n" + s
    stm = []
    for s in (x for x in joined.split("\n") if x.strip()):
        s = s.replace("eft_cache%", "").replace("eft_par_cache%", "").replace("_dl", "")
        s = re.sub(r"(\d)\.(?!\d)", r"\1.0", s)
        stm.append(s)
    return stm


STM = transcribe()


def mu2(env, mpm):
    ns = dict(env)
    ns["sqrt"] = mp.sqrt if mpm else math.sqrt
    for s in STM:
        exec(s, {}, ns)
    return ns["EFT_mu1"], ns["EFT_mu2"]


ZERO = ["EFTOmegaV", "EFTOmegaP", "EFTOmegaPP", "EFTOmegaPPP", "EFTc", "EFTcdot", "EFTcdotdot", "EFTGamma1P",
        "EFTGamma1PP", "EFTGamma2V", "EFTGamma2P", "EFTGamma2PP", "EFTGamma3V", "EFTGamma3P", "EFTGamma3PP",
        "EFTGamma4V", "EFTGamma4P", "EFTGamma4PP"]


def env_from(bg, g, mpm):
    cv = (lambda x: mp.mpf(x)) if mpm else float
    e = {k: cv(0) for k in ZERO}
    e.update({k: cv(v) for k, v in bg.items()})
    e["EFTGamma1V"] = cv(g)
    return e


def main():
    d = np.load(NPZ)
    a_grid = d["eft_a"]
    i = int(np.argmin(np.abs(a_grid - 0.0101)))
    a = float(a_grid[i])
    bg_camb = {k: float(d["eft_" + k][i]) for k in ["adotoa", "Hdot", "Hdotdot", "grhob_t", "grhoc_t"]}
    bg_camb["a"] = a
    import json
    meta = json.loads(str(d["meta"]))
    h0 = meta.get("derived", {}).get("H0") or 67.0
    h0_Mpc = None
    # h0_Mpc = H0/c in 1/Mpc
    h0_Mpc = float(meta.get("eft_params", {}).get("h0_Mpc", 0)) or h0 / 299792.458
    bg_camb["h0_Mpc"] = h0_Mpc
    print(f"a = {a:.6g}, h0_Mpc = {h0_Mpc:.6g}; CAMB mu2/H^2 stored at this a = "
          f"{d['eft_EFT_mu2'][i] * (a / bg_camb['adotoa'])**2:.4g}")
    # consistency check of CAMB background vs matter-only relation:
    H = bg_camb["adotoa"]
    rhom = bg_camb["grhob_t"] + bg_camb["grhoc_t"]
    print(f"(adotoa^2 - Hdot) / (grhom_t/2) = {(H**2 - bg_camb['Hdot']) / (rhom / 2):.6f}  (1 = matter-only)")

    # exact matter + Lambda background with same a, same grhom_t, Lambda part chosen arbitrary but consistent
    x = mp.mpf(rhom)            # a^2 rho_m  (∝ 1/a)
    y = mp.mpf(3) * mp.mpf(H)**2 - x  # a^2 rho_rest treated as Lambda (∝ a^2)
    Hx = mp.sqrt((x + y) / 3)
    Hd = (2 * y - x) / 6
    Hdd = Hx * (4 * y + x) / 6
    bg_exact = {"a": mp.mpf(a), "adotoa": Hx, "Hdot": Hd, "Hdotdot": Hdd, "grhob_t": x, "grhoc_t": mp.mpf(0),
                "h0_Mpc": mp.mpf(h0_Mpc)}

    print(f"{'Gamma1':>8s} {'CAMB-bg dbl':>14s} {'CAMB-bg mp60':>14s} {'exact-bg dbl':>14s} {'exact-bg mp60':>14s}"
          "   (mu_2/H^2)")
    for g in [1.0, 1e-1, 1e-2, 1e-3, 1e-5, 1e-7, 1e-9, 1e-12]:
        row = []
        for bg in (bg_camb, bg_exact):
            for mpm in (False, True):
                try:
                    _, m2 = mu2(env_from(bg, g, mpm), mpm)
                    row.append(float(m2 * (bg["a"] / bg["adotoa"])**2))
                except (ZeroDivisionError, ValueError) as ex:
                    row.append(float("nan"))
        print(f"{g:8.0e} " + " ".join(f"{v:14.5g}" for v in row))
    print("threshold: mu/H^2 < -(rate=10)^2 = -100 -> rejected")


if __name__ == "__main__":
    sys.exit(main())

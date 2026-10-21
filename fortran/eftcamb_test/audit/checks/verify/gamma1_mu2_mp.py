"""Verifier: evaluate the Fortran EFT_mu2 expression (06_abstract_EFTCAMB_model.f90) verbatim in mpmath at high precision
for Pure EFT constant Gamma1 on an LCDM background, with and without radiation, to separate rounding from an analytic
non-cancellation. Units: h0_Mpc = 1, conformal time; grho_x_t = 3 Omega_x a^2 * (a-scaling)."""

import re
from pathlib import Path

import mpmath as mp

mp.mp.dps = 60
SRC = Path(__file__).resolve().parents[4] / "eftcamb" / "06_abstract_EFTCAMB_model.f90"
lines = SRC.read_text().splitlines()
# body of compute_mass_factors from "rho_d =" (second occurrence region, after line 764) to end of EFT_mu2
start = next(i for i in range(760, 800) if lines[i].strip().startswith("rho_d ="))
end = next(i for i in range(873, 910) if lines[i].rstrip().endswith("/(16._dl*F1**2._dl*F2**4._dl*F3**2._dl)") and i > 873)
code = "\n".join(lines[start : end + 1])
code = re.sub(r"!.*", "", code)
code = re.sub(r"&\s*\n\s*", " ", code)
code = code.replace("_dl", "").replace("eft_cache%", "").replace("eft_par_cache%", "")
code = re.sub(r"(\d)\.(?!\d)", r"\1.0", code)
code = code.replace("sqrt(", "mp.sqrt(")
stmts = [s for s in code.split("\n") if s.strip() and not s.strip().startswith("call")]


def mu2(a, G1, Om=0.3, Orad=0.0, h0=1.0):
    a = mp.mpf(a)
    OL = 1 - Om - Orad
    H2 = h0**2 * (Om / a + Orad / a**2 + OL * a**2)  # conformal adotoa^2
    adotoa = mp.sqrt(H2)
    Hdot = h0**2 * (-Om / (2 * a) - Orad / a**2 + OL * a**2)
    dHdot_da = h0**2 * (Om / (2 * a**2) + 2 * Orad / a**3 + 2 * OL * a)
    Hdotdot = adotoa * a * dHdot_da
    env = dict(a=a, adotoa=adotoa, Hdot=Hdot, Hdotdot=Hdotdot, h0_Mpc=mp.mpf(h0), mp=mp,
               grhob_t=mp.mpf(0), grhoc_t=3 * Om * h0**2 / a, EFTc=0, EFTcdot=0, EFTcdotdot=0,
               EFTOmegaV=0, EFTOmegaP=0, EFTOmegaPP=0, EFTOmegaPPP=0,
               EFTGamma1V=mp.mpf(G1), EFTGamma1P=0, EFTGamma1PP=0,
               EFTGamma2V=0, EFTGamma2P=0, EFTGamma2PP=0,
               EFTGamma3V=0, EFTGamma3P=0, EFTGamma3PP=0, EFTGamma4V=0, EFTGamma4P=0, EFTGamma4PP=0)
    for s in stmts:
        exec(s.strip(), {"mp": mp}, env)
    return env["EFT_mu2"] * (a / adotoa) ** 2, env


if __name__ == "__main__":
    for orad in [0.0, 1e-4]:
        for dps in [60, 16]:
            mp.mp.dps = dps
            print(f"Orad={orad} dps={dps}")
            for G1 in [1e-2, 1e-3, 1e-5, 1e-7]:
                m, env = mu2(0.01, G1, Orad=orad)
                print(f"  G1={G1:8.0e} mu2/H^2={mp.nstr(m, 6):>14s}  3F2^2+F1F3={mp.nstr(3*env['F2']**2+env['F1']*env['F3'], 6)}")

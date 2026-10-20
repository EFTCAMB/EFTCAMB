"""(iv) sanity scale of the mass eigenvalues for canonical quintessence (5_quint_*): the decoupling-limit pi mass of a
canonical scalar (delta phi = phidot pi) is m_pi^2 = -3 dH/dt - 3 H phiddot/phidot = -3 dH/dt - (3/2) H cdot/c,
i.e. m_pi^2/H^2 = -3 (H'/H^2 - 1) - (3/2) EFTcdot/(H EFTc) in EFTCAMB conformal variables (EFTc = c a^2).
Compare with the code's EFT_mu1, EFT_mu2 (in units of H^2 = (adotoa/a)^2, as in the mass check)."""
import sys
import numpy as np
sys.path.insert(0, ".")
import t6lib, harness
g = np.logspace(-2, 0, 201)
for m in ("5_quint_exponential", "5_quint_inverse_power", "4_Des_mc5e_wCDM_1"):
    ini = harness.PARAMS_DIR / f"{m}.ini"
    r = t6lib.probe("A0", ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0",
                                          "EFT_mass_stability": "T"}, a_grid=g, feedback=0)
    a, H = r["a"], r["adotoa"]
    mpi = -3 * (r["Hdot"] / H**2 - 1) - 1.5 * r["EFTcdot"] / (H * r["EFTc"])
    mu1 = r["EFT_mu1"] * (a / H) ** 2
    mu2 = r["EFT_mu2"] * (a / H) ** 2
    q = lambda x: f"[{np.nanmin(x): .3e}, {np.nanmax(x): .3e}] median {np.nanmedian(x): .3e}"
    print(f"{m}: decoupling m_pi^2/H^2 in {q(mpi)}")
    print(f"{'':{len(m)}s}  code mu1/H^2 in {q(mu1)}")
    print(f"{'':{len(m)}s}  code mu2/H^2 in {q(mu2)}   (mass check rejects if < -rate^2 = -100)")

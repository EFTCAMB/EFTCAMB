"""Evaluate the code's mu1/mu2 formula on a dust-only-consistent background: Hdot, Hdotdot replaced by
H^2 - c - a^2 rho_d/2 and its derivative (rho_d = b+c dust, as the formula assumes). Separates the early-time
1e9 values (background inconsistency: radiation + massive nu not in rho_d, divided by powers of c) from the
late-time O(-3e3) values."""
import sys
import numpy as np
import mpmath
import massverdict_indep as M
mpmath.mp.dps = 50
M.WORKER = M.WORKER.replace("grid = np.array([0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0])", "grid = np.logspace(-2, 0, 9)")
stmts = M.translate()
for model in sys.argv[2:]:
    r = M.dump(sys.argv[1], model)
    print(model)
    for c in r["rows"]:
        a, H = c["a"], c["adotoa"]
        cc = dict(c)
        dust = 0.5 * (c["grhob_t"] + c["grhoc_t"])  # = a^2 rho_d / 2
        cc["Hdot"] = H**2 - c["EFTc"] - dust
        cc["Hdotdot"] = 2 * H * cc["Hdot"] - c["EFTcdot"] + H * dust  # d(a^2 rho_d)/dtau = -H a^2 rho_d
        ns = M.evaluate(stmts, cc, r["h0_Mpc"], mp=mpmath)
        n = (a / H) ** 2
        print(f"  a={a:.4f} c/H^2={c['EFTc']/H**2:.3g}  code mu1,mu2/H^2 = {c['EFT_mu1']*n:.4g}, {c['EFT_mu2']*n:.4g}"
              f"   dust-consistent mu1,mu2/H^2 = {float(ns['C_EFT_mu1'])*n:.4g}, {float(ns['C_EFT_mu2'])*n:.4g}")

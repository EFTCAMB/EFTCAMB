"""Dense-grid companion of massverdict_indep.py: where is min mu/H^2, and does it survive if rho_d includes massive nu?"""
import sys, json
import numpy as np
import massverdict_indep as M
import mpmath
mpmath.mp.dps = 40
M.WORKER = M.WORKER.replace("grid = np.array([0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0])", "grid = np.logspace(-2, 0, 201)")
tree = sys.argv[1]
stmts = M.translate()
for model in sys.argv[2:]:
    r = M.dump(tree, model)
    out = []
    for c in r["rows"]:
        a, H = c["a"], c["adotoa"]
        if a <= 1e-8:
            continue
        norm = (a / H) ** 2
        full = 0.5 * (c["grhob_t"] + c["grhoc_t"]) + (2 / 3) * (c["grhog_t"] + c["grhor_t"]) + 0.5 * (c["grhonu_tot"] + c["gpinu_tot"])
        rf = 2 * full / a ** 2
        try:
            mpf = M.evaluate(stmts, c, r["h0_Mpc"], mp=mpmath, override={"rho_d": lambda ns, _r=rf: mpmath.mpf(_r)})
            m1f, m2f = float(mpf["C_EFT_mu1"]) * norm, float(mpf["C_EFT_mu2"]) * norm
        except ZeroDivisionError:
            m1f = m2f = float("nan")
        out.append((a, c["EFTc"] / H**2, c["EFT_mu1"] * norm, c["EFT_mu2"] * norm, m1f, m2f))
    o = np.array(out)
    mn = np.fmin(o[:, 2], o[:, 3]); mnf = np.fmin(o[:, 4], o[:, 5])
    i = np.nanargmin(mn); j = np.nanargmin(mnf)
    print(f"{model}: code min mu/H^2 = {mn[i]:.4g} at a={o[i,0]:.4f} (c/H^2={o[i,1]:.3g}); rejected points (< -100): "
          f"{np.sum(mn < -100)}/{len(mn)}, NaN {np.sum(np.isnan(mn))}")
    print(f"   rho_d incl. nu+rad: min mu/H^2 = {mnf[j]:.4g} at a={o[j,0]:.4f}; rejected points: {np.sum(mnf < -100)}/{len(mnf)}")
    for k in (0, 50, 100, 150, 200):
        if k < len(o):
            print("   a=%.4f c/H2=%.3g mu1=%.4g mu2=%.4g | full-rho mu1=%.4g mu2=%.4g" % tuple(o[k]))

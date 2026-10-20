"""(iv) EFT_mass_stability: min(mu1, mu2)/H^2 after turn-on on A0 vs the candidate derivative fix (fixaddl), and the
number of points below -rate^2 (rate = EFT_mass_stability_rate = 10), for canonical / healthy test models."""
import sys
import numpy as np
sys.path.insert(0, ".")
import t6lib, harness
FIX = str(t6lib.WORK / "fixaddl")
models = sys.argv[1:] or ["4_Des_mc5e_wCDM_1", "5_quint_exponential", "4_Des_F_R_LCDM_1", "2_PEFT_Omega_const_1",
                          "3_RPH_MK_const_1", "5_Kmouflage_1", "2_PEFT_wDE_wCDM_1"]
g = np.logspace(-2, 0, 401)
for m in models:
    ini = harness.PARAMS_DIR / f"{m}.ini"
    for tree in ("A0", FIX):
        r = t6lib.probe(tree, ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0",
                                              "EFT_mass_stability": "T"}, a_grid=g, feedback=0)
        a, H = r["a"], r["adotoa"]
        sel = a > r["pert_turn_on"]
        x1 = r["EFT_mu1"][sel] * (a[sel] / H[sel]) ** 2
        x2 = r["EFT_mu2"][sel] * (a[sel] / H[sel]) ** 2
        x = np.minimum(x1, x2)
        nnan = int((~np.isfinite(x1) | ~np.isfinite(x2)).sum())
        fin = np.isfinite(x)
        nbad = int((x[fin] < -r["EFT_mass_stability_rate"] ** 2).sum())
        i = np.nanargmin(np.where(fin, x, np.nan)) if fin.any() else None
        print(f"{m:24s} {('A0' if tree == 'A0' else 'fixaddl'):8s} turn_on={r['pert_turn_on']:.0e} "
              f"min mu/H^2={x[i] if i is not None else float('nan'): .3e} at a={a[sel][i] if i is not None else float('nan'):.3f} "
              f"median mu1/H^2={np.nanmedian(x1): .3e} mu2/H^2={np.nanmedian(x2): .3e} "
              f"pts<-rate^2: {nbad}/{int(sel.sum())} NaN: {nnan}")

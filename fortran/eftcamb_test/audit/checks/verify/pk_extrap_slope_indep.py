"""Independent check of pk-extrap-slope-bias: same A0 data, compare slope formulas and the k grid."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import harness
T = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0"
sys.path.insert(0, T)
import numpy as np
with harness.slot(1):
    import camb
    print("camb from", camb.__file__)
    def pars(kmax, nl, model=None):
        p = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2.1e-9, ns=0.965, WantCls=False,
                            WantTransfer=True, kmax=kmax, redshifts=[0.0], NonLinear=nl)
        return p
    for nl in ("NonLinear_none", "NonLinear_both"):
        isnl = nl != "NonLinear_none"
        lo = camb.get_results(pars(2.0, nl)); hi = camb.get_results(pars(100.0, nl))
        kh, z, pk = lo.get_linear_matter_power_spectrum(nonlinear=isnl)
        lk = np.log(kh); lp = np.log(pk[0])
        dl = np.diff(lk)
        print(nl, "nk", len(kh), "kh_max", kh[-1], "last 10 dlogk", np.round(dl[-10:], 5))
        s_end = (lp[-1] - lp[-2]) / dl[-1]
        s_mean = np.mean(np.diff(lp[-10:])) / dl[-1]
        s_secant = (lp[-1] - lp[-10]) / (lk[-1] - lk[-10])
        print("  slopes: end %.4f  EFT-mean %.4f  true secant(last9) %.4f" % (s_end, s_mean, s_secant))
        Ptrue = hi.get_matter_power_interpolator(nonlinear=isnl, silent=True)
        r = np.array([1.5, 2, 5, 10]); k = kh[-1] * r
        pt = Ptrue.P(0.0, k)
        for name, s in (("end", s_end), ("mean", s_mean)):
            pe = pk[0][-1] * (r ** s)
            print("  pure power law %-4s rel err" % name, np.round(pe / pt - 1, 3))
        pi = lo.get_matter_power_interpolator(nonlinear=isnl, extrap_kmax=30 * kh[-1], silent=True).P(0.0, k)
        print("  A0 interpolator rel err      ", np.round(pi / pt - 1, 3))
        # local true slope of hi P at kh_max and beyond
        lh = Ptrue
        for kk in (kh[-1], 2 * kh[-1], 10 * kh[-1]):
            e = 1e-3
            print("  true local slope at k=%.2f: %.3f" % (kk, (np.log(lh.P(0, kk * (1 + e))) - np.log(lh.P(0, kk * (1 - e)))) / (2 * e)))

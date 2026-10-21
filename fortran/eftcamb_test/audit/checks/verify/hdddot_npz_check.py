"""Independent numeric check on the harness reference npz (A0 and B1): Fortran cache Hdotdotdot vs (i) literal
python transcription of 06p1:141-168, (ii) proposed fix, (iii) finite-difference d(Hdotdot)/dtau (np.gradient)."""
import sys
import numpy as np
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks/verify")
from hdddot_sympy import fortran, proposed  # noqa: E402
R = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/runs"
for tree in ("A0", "B1"):
    for m in ("5_Kmouflage_1", "5_Kmimic_1", "5_ADE_1", "5_quint_exponential"):
        try:
            d = np.load(f"{R}/{tree}/{m}.npz", allow_pickle=True)
        except FileNotFoundError:
            print(tree, m, "missing"); continue
        a = d["eft_a_grid"]; g = lambda k: d["eft_" + k]
        args = [a, g("adotoa"), g("Hdot"), g("Hdotdot"), g("gpresm_t"), g("gpresdotm_t"), g("gpresdotdotm_t"),
                g("EFTLambda"), g("EFTLambdadot"), g("EFTLambdadotdot"), g("EFTOmegaV"), g("EFTOmegaP"),
                g("EFTOmegaPP"), g("EFTOmegaPPP"), g("EFTOmegaPPPP")]
        Hddd = g("Hdotdotdot")
        lit, fix = fortran(*args), proposed(*args)
        num = g("adotoa") * np.gradient(g("Hdotdot"), np.log(a))
        sel = (a > 1e-3) & (a < 0.95)
        s = np.max(np.abs(num[sel]))
        r = lambda x, y: np.max(np.abs(x[sel] - y[sel])) / s
        print(f"{tree} {m:20s} maxOmegaP={np.max(np.abs(g('EFTOmegaP')[sel])):.2e} cache-vs-literal {r(Hddd, lit):.1e} "
              f"cache-vs-FD {r(Hddd, num):.2e} fix-vs-FD {r(fix, num):.2e}")

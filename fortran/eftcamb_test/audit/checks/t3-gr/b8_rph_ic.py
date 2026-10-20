"""t3-gr (b): mechanism of the single-k P(k) spikes of RPH near GR (b5): evaluate the pi-equation coefficients and pi
at the perturbation turn-on (a = 0.01) for k on and around the spike (RPHmassP0 = 1e-6, spike at k = 0.1726/Mpc).
Hypothesis: the standard EFT IC (09_EFTCAMB_IC.f90 case 1) pi_ini = -E/(C + k^2 D1 + k^4 D2) is evaluated where the
denominator nearly vanishes.
"""

import numpy as np

from common import PARAMS, RUNS, harness

INI = PARAMS / "3_RPH_MK_const_1.ini"
OUT = RUNS / "b8" / "RPH_m1e-6_evol.npz"
ref = harness.load(RUNS / "b1" / "REF_base__A0g.npz")
kk = ref["pk_k"]
i0 = int(np.argmin(abs(kk - 0.1726)))
KS = [float(k) for k in kk[i0 - 3:i0 + 4]]
if not OUT.exists():
    harness.run_model("A0", INI, OUT, overrides={"RPHmassP0": "1e-06", "get_tensor_cls": "F", "do_lensing": "F",
                                                "l_max_scalar": "100"}, evol_k=KS, threads=1, timeout=1500)
d = harness.load(OUT)
a = np.atleast_2d(d["evol_a"])[0]
j = int(np.searchsorted(a, 0.01))  # first sampled a after the turn-on
print(f"a = {a[j]:.4g}")
full = harness.load(RUNS / "b1" / "RPH__s1em06.npz")  # same model, full settings, same k grid as ref
r = np.abs(full["pk"][0] / ref["pk"][0] - 1)
for i, k in enumerate(np.atleast_2d(d["evol_k"])[:, 0]):
    C, D1, D2, E = (d[f"evol_{f}"][i, j] for f in ("EFTpiC", "EFTpiD1", "EFTpiD2", "EFTpiE"))
    den = C + k * k * D1 + k ** 4 * D2
    pi = d["evol_pi"][i, j]
    print(f"k={k:.5f} C={C:+.4e} k2D1={k * k * D1:+.4e} den={den:+.3e} E={E:+.3e} -E/den={-E / den if den else 0:+.3e} "
          f"pi={pi:+.3e} " + (f"|dP/P|={r[i0 - 3 + i]:.1e}" if r is not None else ""))

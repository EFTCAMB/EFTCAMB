"""Is the 'x626 at turn-on 1e-6' statement evidence against IC2?  Compare IC1 and IC2 (turn-on 1e-6) with the turn-on
1e-8 runs (where IC1 and IC2 agree) for 2_PEFT_Omega_const_2.  Inputs: the npz of ic23_early.py."""

import sys

import numpy as np

sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

D = harness.AUDIT_ROOT / "work" / "verify" / "ic23" / (sys.argv[1] if len(sys.argv) > 1 else "early")
r = {n: harness.load(D / f"{n}.npz") for n in ("ic1_ton1e-8", "ic2_ton1e-8", "ic1_ton1e-6", "ic2_ton1e-6")}
k = r["ic1_ton1e-8"]["pk_k"]
ref = r["ic1_ton1e-8"]["pk"][0]
q = r["ic2_ton1e-6"]["pk"][0] / r["ic1_ton1e-6"]["pk"][0]
i626 = int(np.argmax(np.abs(q - 1)))
print(f"max |P_IC2/P_IC1 - 1| at turn-on 1e-6 = {abs(q[i626]-1):.4g} at k={k[i626]:.4g}")
print(f"{'k':>8} {'IC1(1e-6)/ref':>14} {'IC2(1e-6)/ref':>14} {'IC2(1e-8)/ref':>14} {'IC2/IC1(1e-6)':>14}")
for kk in (0.1, 0.5, 1, 2, 5, 8, 10, 15, 20, 30, k[i626]):
    i = int(np.argmin(abs(k - kk)))
    print(f"{k[i]:8.4g} {r['ic1_ton1e-6']['pk'][0][i]/ref[i]:14.4g} {r['ic2_ton1e-6']['pk'][0][i]/ref[i]:14.4g} "
          f"{r['ic2_ton1e-8']['pk'][0][i]/ref[i]:14.4g} {q[i]:14.4g}")
for n in ("ic1_ton1e-6", "ic2_ton1e-6", "ic2_ton1e-8"):
    c = harness.compare(r[n], r["ic1_ton1e-8"])
    print(n, "vs ic1_ton1e-8: TT %.2e EE %.2e TE %.2e  max|dP/P| %.3e  sigma8 %.6f vs %.6f" % (
        c["cl"]["total"]["TT"], c["cl"]["total"]["EE"], c["cl"]["total"]["TE"], c["pk"][0],
        float(r[n]["sigma8"][0]), float(r["ic1_ton1e-8"]["sigma8"][0])))

"""Which pi IC is right in the radiation era?  Evolve Omega_const_2 (Omega const = 0.1) from turn-on 1e-8 and read pi at
a = 1e-6, compare with the IC-2 series formula and the IC-1 attractor at a = 1e-6 (from the IC23DUMP of the 1e-6 run).
Also compare P(k) of the 1e-6-turn-on runs (IC1, IC2) with the 1e-8-turn-on runs."""

import json
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "verify" / "ic23"
SW = W / "A0_sw"
INI = harness.PARAMS_DIR / "2_PEFT_Omega_const_2.ini"
OUT = W / "early"
OUT.mkdir(exist_ok=True)
KS = [0.01, 0.1, 0.5436, 4.255, 10.0]
runs = {f"ic{ic}_ton{t}": ({"EFT_IC_type": str(ic), "EFTCAMB_turn_on_time": t}, {"IC23_DUMP": "1"})
        for ic in (1, 2) for t in ("1e-8", "1e-6")}


def one(item):
    name, (ov, env) = item
    return name, harness.run_model(SW, INI, OUT / f"{name}.npz", overrides=ov, extra_env=env, evol_k=KS, threads=1)


with ThreadPoolExecutor(4) as ex:
    st = dict(ex.map(one, runs.items()))
print({k: v["status"] for k, v in st.items()})
d = {k: harness.load(OUT / f"{k}.npz") for k in runs}
a = d["ic1_ton1e-8"]["evol_a"]
print("evol a[0] =", a[0])
h0 = 70.0 / 299792.458  # placeholder, replaced below by H0 from the run
H0 = d["ic1_ton1e-8"]["meta"]["derived"].get("H0", None)
for name in runs:
    pi = d[name]["evol_pi"]
    tau = d[name]["evol_tau"][0, 0] if "evol_tau" in d[name] else None
    print(name, "pi(a=1e-6) per k:", " ".join(f"{x: .4e}" for x in pi[:, 0]), " tau(1e-6)=", tau)
# P(k) comparisons
k = d["ic1_ton1e-8"]["pk_k"]
for ref in ("ic1_ton1e-8", "ic2_ton1e-8"):
    for name in runs:
        if name == ref:
            continue
        r = d[name]["pk"][0] / d[ref]["pk"][0] - 1
        cl = harness.compare(d[name], d[ref])["cl"]["total"]
        sel = [np.argmin(abs(k - kk)) for kk in (0.1, 1, 5, 10, 20)]
        print(f"{name:>13} vs {ref}: max|dP/P|={np.max(abs(r)):.3e} at k={k[np.argmax(abs(r))]:.3g}; "
              + " ".join(f"k={k[i]:.3g}:{r[i]:.2e}" for i in sel) + f"  TT {cl['TT']:.2e} EE {cl['EE']:.2e}")

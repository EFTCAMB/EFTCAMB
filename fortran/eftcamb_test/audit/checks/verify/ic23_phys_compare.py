import sys, json
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness, numpy as np
OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic23phys"
model = sys.argv[1]; ref = sys.argv[2]; others = sys.argv[3:]
R = harness.load(str(OUT / f"{model}_{ref}.npz"))
def info(D):
    m = D["meta"] if isinstance(D["meta"], dict) else json.loads(str(D["meta"]))
    p = m.get("eft_params", {})
    return {k: v for k, v in p.items() if "IC" in k or "turn" in k}
print("ref", ref, info(R), "sigma8", R["sigma8"])
for o in others:
    D = harness.load(str(OUT / f"{model}_{o}.npz"))
    c = harness.compare(R, D)
    a = D["evol_a"][0]; i = int(np.argmax(a > 0.0102))
    pi = D["evol_pi"][:, i] if "evol_pi" in D else None
    pir = R["evol_pi"][:, i] if "evol_pi" in R else None
    print(o, info(D), "sigma8", D["sigma8"], "worst_cl", c["worst_cl"], "worst_pk", c["worst_pk"])
    print("   a=%.4f pi(ref)=%s pi(this)=%s" % (a[i], np.array2string(pir, precision=3), np.array2string(pi, precision=3)))
    k = D["pk_k"]; r = D["pk"][0] / R["pk"][0] - 1
    for kk in (0.01, 0.1, 1.0, 5.0):
        j = np.argmin(abs(k - kk)); print("   dP/P(k=%g, z=%g) = %.3e" % (k[j], D["pk_z"][0], r[j]))

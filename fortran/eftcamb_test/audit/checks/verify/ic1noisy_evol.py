"""Verifier: inspect early-time pi coefficients and pi evolution for 5_Kmimic_1 / 5_hdsk_TG_2 on a fine a grid."""
import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
REPO = "/Users/marco/IW/eftcamb_update/EFTCAMB"
sys.path.insert(0, REPO + "/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402
import numpy as np  # noqa: E402

tree, model, turnon, env = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
for e in env:
    k, v = e.split("=")
    os.environ[k] = v
sys.path.insert(0, str(harness.tree_path(tree)))
import camb  # noqa: E402

out = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy"
ini = out / f"{model}_{turnon}.ini"
ini.write_text(f"DEFAULT({harness.PARAMS_DIR}/{model}.ini)\nEFTCAMB_turn_on_time = {turnon}\nEFTCAMB_write_background=F\n"
               "")
with harness.slot(1):
    cp = camb.read_ini(str(ini), no_validate=True)

    res = camb.get_results(cp); print("results ok", flush=True)
    t = float(turnon)
    a = np.logspace(np.log10(t) - 0.02, np.log10(t) + 2, 41)
    ks = [1.0, 5.0]
    f, arr = cp.EFTCAMB.get_scale_evolution(res, ks, a)
names = arr.dtype.names
print([n for n in names if "pi" in n.lower()][:40])
for ik, k in enumerate(ks):
    print("k =", k)
    print("%10s %11s %11s %11s %11s %11s %11s %11s" % ("a", "C+k2D1+k4D2", "E", "piQS", "pi", "pidot", "G1", "A1"))
    for ia in range(0, len(a), 2):
        r = arr[ik, ia]
        den = r["EFTpiC"] + k**2 * r["EFTpiD1"] + k**4 * r["EFTpiD2"]
        qs = -cp.H0 / 2.99792458e5 * r["EFTpiE"] / den if den != 0 else np.nan
        print("%10.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e" % (
            a[ia], den, r["EFTpiE"], qs, r["pi"], r["pidot"], r["EFTGamma1V"], r["EFTpiA1"]))

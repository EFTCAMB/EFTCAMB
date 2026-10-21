"""Independent check of get-evolution-pi-highk-tolerance: convergence ladder in IntTolBoost, lAccuracyBoost effect,
upstream get_time_evolution path (Weyl, delta_cdm) for EFT model and GR."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "2"
TREE = sys.argv[1]
INI = sys.argv[2]
EFTFLAG0 = len(sys.argv) > 3 and sys.argv[3] == "gr"
sys.path.insert(0, TREE)
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import tempfile
import numpy as np
import harness
import camb
tmp = tempfile.mkdtemp()
fast = os.path.join(tmp, "model.ini")
open(fast, "w").write(f"DEFAULT({INI})\nl_max_scalar = 1500\nget_tensor_cls = F\nget_transfer = F\n"
                      "EFTCAMB_write_background = F\nfeedback_level = 0\n" + ("EFTflag = 0\n" if EFTFLAG0 else ""))
k = 1.3172836143
a = np.logspace(np.log10(0.02), 0, 400)
res = {}
with harness.slot(2):
    for tag, itb, lab in [("d", 1, 1), ("t2", 2, 1), ("t4", 4, 1), ("t8", 8, 1), ("t6", 6, 1), ("l4", 1, 4)]:
        p = camb.read_ini(fast)
        p.Accuracy.IntTolBoost = itb
        p.Accuracy.lAccuracyBoost = lab
        r = camb.get_results(p)
        z = 1 / a - 1
        eta = r.conformal_time(z)
        tev = r.get_time_evolution(k, eta, ["Weyl", "delta_cdm", "delta_baryon"], lAccuracyBoost=1)
        tevc = r.get_time_evolution(k, eta[::4], ["Weyl", "delta_cdm"], lAccuracyBoost=1)
        d = {"weyl": tev[:, 0], "dc": tev[:, 1], "weyl_c": tevc[:, 0], "dc_c": tevc[:, 1]}
        if not EFTFLAG0:
            f = p.EFTCAMB.get_scale_evolution(r, k, a)[1][0]
            fc = p.EFTCAMB.get_scale_evolution(r, k, a[::4])[1][0]
            for n in ("pi", "Psi", "Phi"):
                if n in f.dtype.names:
                    d[n] = f[n].copy(); d[n + "_c"] = fc[n].copy()
        res[tag] = d
ref = res["t8"]
print("model:", os.path.basename(INI), "GR" if EFTFLAG0 else "EFT")
for name in [n for n in res["d"] if not n.endswith("_c")]:
    m = np.max(np.abs(ref[name]))
    row = []
    for tag in ("d", "t2", "t4", "t6", "l4"):
        row.append(f"{tag}:{np.max(np.abs(res[tag][name] - ref[name])) / m:.2e}")
    sc = np.max(np.abs(res["d"][name][::4] - res["d"][name + "_c"])) / m
    sct = np.max(np.abs(res["t8"][name][::4] - res["t8"][name + "_c"])) / m
    print(f"{name:6s} vs IntTolBoost8: " + " ".join(row) + f" | stops-dependence default {sc:.2e}, t8 {sct:.2e}")

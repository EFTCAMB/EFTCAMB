"""Independent check of ic1-pi-attractor-pole: dense k scan around k_c, test pi_ini*(k-k_c) ~ const (simple pole),
pi_ini vs -h0*E/den, and B1/A0 reference P(k) near k_c."""
import json, os, subprocess, sys
from pathlib import Path
import numpy as np
REPO = Path("/Users/marco/IW/eftcamb_update/EFTCAMB")
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness
OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1pole"
OUT.mkdir(parents=True, exist_ok=True)
MODEL = "3_RPH_MB_const_1"
EVOL = r'''
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_results(cp)
ks = json.loads(sys.argv[3])
a = np.array([0.0100001, 0.01001, 1.0])
f, ev = cp.EFTCAMB.get_scale_evolution(res, ks, a)
out = {n: [[float(ev[n][i][j]) for j in range(3)] for i in range(len(ks))] for n in ("EFTpiC","EFTpiD1","EFTpiD2","EFTpiE","pi","clxc")}
out["h0"] = float(cp.H0/2.99792458e5)
print(json.dumps(out))
'''
ini = OUT / "w.ini"
ini.write_text(f"DEFAULT({harness.PARAMS_DIR / (MODEL + '.ini')})\noutput_root = {OUT / 'x'}\n"
               "EFTCAMB_write_background = F\nget_tensor_cls = F\nnum_redshiftwindows = 0\ndo_lensing = F\n"
               "l_max_scalar = 300\nget_transfer = F\nfeedback_level = 0\n")
(OUT / "evol.py").write_text(EVOL)
ks = list(np.concatenate([np.linspace(1.55e-3, 1.70e-3, 31), [1.0e-3, 3.0e-3]]))
with harness.slot(1):
    r = subprocess.run([sys.executable, str(OUT / "evol.py"), str(harness.tree_path("A0")), str(ini), json.dumps(ks)],
                       cwd=OUT, env=dict(os.environ, OMP_NUM_THREADS="1"), capture_output=True, text=True, timeout=600)
d = json.loads(r.stdout.strip().splitlines()[-1])
h0 = d["h0"]
C = np.array(d["EFTpiC"])[:, 0]; D1 = np.array(d["EFTpiD1"])[:, 0]; D2 = np.array(d["EFTpiD2"])[:, 0]
E = np.array(d["EFTpiE"])[:, 0]; pi = np.array(d["pi"])[:, 0]; pi2 = np.array(d["pi"])[:, 1]; dc = np.array(d["clxc"])[:, 2]
k = np.array(ks)
den = C + k**2 * D1 + k**4 * D2
kc = np.sqrt(-C[0] / D1[0])
print(f"h0={h0:.4e} kc={kc:.6e}  C={C[0]:.3e} D1={D1[0]:.3e} D2={D2[0]:.3e}")
for i in range(len(k)):
    print(f"k={k[i]:.6f} den={den[i]:+.3e} pi(a=1.00001e-2)={pi[i]:+.3e} -h0E/den={-h0*E[i]/den[i]:+.3e} "
          f"pi*(k-kc)={pi[i]*(k[i]-kc):+.3e} pi(a=.01001)={pi2[i]:+.3e} dc(1)={dc[i]:+.3e}")
for tr in ("A0", "B1"):
    p = harness.AUDIT_ROOT / "runs" / tr / f"{MODEL}.npz"
    g = harness.AUDIT_ROOT / "runs" / tr / "1_EFT_GR.npz"
    if p.exists() and g.exists():
        a, b = harness.load(p), harness.load(g)
        kk = a["pk_k"]; sel = (kk > 1.2e-3) & (kk < 2.6e-3)
        pg = np.exp(np.interp(np.log(kk[sel]), np.log(b["pk_k"]), np.log(b["pk"][0])))
        print(tr, " ".join(f"{x:.5f}:{y:.3g}" for x, y in zip(kk[sel], a["pk"][0][sel] / pg)))

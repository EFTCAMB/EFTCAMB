"""Independent check of QS_sigma vs analytic f(R) quasi-static Sigma = 1/(1+Omega) (Omega = f_R).
Uses harness.run_model (evol_k) on A0, 4_Des_F_R_LCDM_1; compares QS_sigma, QS_mu with the analytic
f(R) expressions Sigma = 1/(1+Omega), mu(k->inf) = 4/3 /(1+Omega)."""
import sys, numpy as np
REPO = "/Users/marco/IW/eftcamb_update/EFTCAMB"
sys.path.insert(0, REPO + "/fortran/eftcamb_test/audit/checks")
import harness
SP = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
out = SP + "/work/verify/qs_sigma/fR1.npz"
import os
if not os.path.exists(out): st = harness.run_model(SP + "/trees/A0", REPO + "/fortran/eftcamb_test/parameters/4_Des_F_R_LCDM_1.ini", out,
                       evol_k=[0.01, 0.1, 1.0, 10.0], threads=1)
pass
d = harness.load(out)
A = d["evol_a"]
om = d["evol_EFTOmegaV"]
for i, k in enumerate(d["evol_k"][:, 0]):
    a = A[i]
    sel = a > 0.2
    sig_an = 1 / (1 + om[i, sel])
    print(f"k={k}: max|QS_sigma-QS_mu|={np.max(np.abs(d['evol_QS_sigma'][i,sel]-d['evol_QS_mu'][i,sel])):.3e}"
          f"  max|QS_sigma-1/(1+Om)|={np.max(np.abs(d['evol_QS_sigma'][i,sel]-sig_an)):.3e}"
          f"  max|QS_sigma_gen-1/(1+Om)|={np.max(np.abs(d['evol_QS_sigma_gen'][i,sel]-sig_an)):.3e}"
          f"  max|QS_mu-4/3/(1+Om)|={np.max(np.abs(d['evol_QS_mu'][i,sel]-4/3*sig_an)):.3e}")
j = np.argmin(abs(A[0] - 1))
print("a=1: Omega=", om[:, j], "QS_sigma=", d["evol_QS_sigma"][:, j], "1/(1+Om)=", 1 / (1 + om[:, j]))

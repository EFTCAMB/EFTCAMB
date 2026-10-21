"""Independent check of minkowski-tachpos-lambdadotdot-zero using the A0 reference npz (get_eft_functions path).
Validates the X-dot convention on Lambda->Lambdadot, then rebuilds tachpos with FD Lambdadotdot."""
import sys
import numpy as np
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
R = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/runs/A0"
for m in ("5_quint_exponential", "5_quint_inverse_power"):
    d = harness.load(f"{R}/{m}.npz")
    a = d["eft_a_grid"]; H = d["eft_adotoa"]; L = d["eft_EFTLambda"]; Ld = d["eft_EFTLambdadot"]
    s = (a > 1e-3)
    la = np.log(a)
    Ld_fd = H * np.gradient(L, la, edge_order=2) - 2 * H * L
    Ldd_fd = H * np.gradient(Ld, la, edge_order=2) - 2 * H * Ld
    rel = np.max(np.abs(Ld_fd - Ld)[s][2:-2]) / np.max(np.abs(Ld[s]))
    D1 = -2 * d["eft_EFTc"]  # gamma_i = 0 reduction of D1 (cache ghostpos not filled on this path)
    tach_code = tach0_c = None or (-(a**2) * (-H * Ld) / D1)  # what the code computes with Ldd=0
    tach0 = -(a**2) * (-H * Ld) / D1
    tach_fd = -(a**2) * (Ldd_fd - H * Ld) / D1
    ss = np.where(s)[0][2:-2]
    print(m, "convention check max|Ld_fd-Ld|/max|Ld| =", f"{rel:.2e}",
          "| max|EFTLambdadotdot| in cache =", f"{np.max(np.abs(d['eft_EFTLambdadotdot'])):.2e}",
          "| max|gamma2,3| =", f"{np.max(np.abs(d['eft_EFTGamma2V']))+np.max(np.abs(d['eft_EFTGamma3V'])):.1e}")
    print("   tachpos code vs Ldd=0 formula max rel diff:",
          f"{np.max(np.abs(tach_code-tach0)[ss]/np.abs(tach0[ss])):.2e}")
    print("   points a>1e-3 failing (-tachpos<0): code", int((-tach_code[ss] < 0).sum()), "FD", int((-tach_fd[ss] < 0).sum()),
          "of", len(ss), "; ghost fails", int((-D1[ss] <= 0).sum()))
    print("   |Ldd_fd|/|H Ld| median:", f"{np.median(np.abs(Ldd_fd[ss])/np.abs(H[ss]*Ld[ss])):.2f}")
    i = np.argmin(abs(a - 0.5)); print(f"   a={a[i]:.3f}: tach code {tach_code[i]:.3e} FD {tach_fd[i]:.3e}")

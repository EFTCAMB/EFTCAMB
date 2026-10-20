"""Independent check: with gamma3=gamma4=0 the pi-equation coefficient EFTpiA1 implies
EFT_kinetic = 36 (1+Omega)^2 A1 exactly (expanding (H Omega' + h0 gamma2)^2 gives a single-h0 cross term).
Uses the A0 reference npz (default settings) of each model."""
import sys
import numpy as np
R = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/runs/A0/"
for m in sys.argv[1:]:
    d = np.load(R + m + ".npz", allow_pickle=True)
    a = d["eft_a_grid"]
    O, g3, g4 = d["eft_EFTOmegaV"], d["eft_EFTGamma3V"], d["eft_EFTGamma4V"]
    A1, K, OP, g2 = d["eft_EFTpiA1"], d["eft_EFT_kinetic"], d["eft_EFTOmegaP"], d["eft_EFTGamma2V"]
    sel = (a > 0.05) & (np.abs(K) > 0)
    r = 36 * (1 + O[sel]) ** 2 * A1[sel] / K[sel]
    print(f"{m}: max|g3|={np.abs(g3).max():.1e} max|g4|={np.abs(g4).max():.1e} max|g2*OP|={np.abs(g2*OP).max():.2e} "
          f"ratio 36(1+O)^2A1/K at a=0.1,0.5,1: "
          + " ".join(f"{np.interp(x, a[sel], r):.4f}" for x in (0.1, 0.5, 1.0)) + f"  min={r.min():.4f} max={r.max():.4f}")

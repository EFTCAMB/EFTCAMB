"""Independent check of fullmap-higher-derivs-unset: from the harness reference npz (A0, default settings),
compare each higher-derivative cache field with a finite difference (np.gradient in ln a) of its base field.
No CAMB runs. Prints max|field|, max|FD|, rel err over a in [1e-3, 1]."""
import glob
import os
import numpy as np

R = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/runs/A0"
models = sorted(os.path.basename(p)[:-4] for p in glob.glob(R + "/[3-5]_*.npz"))
for m in models:
    d = np.load(f"{R}/{m}.npz", allow_pickle=True)
    a = d["eft_a_grid"]; la = np.log(a); H = d["eft_adotoa"]; s = a > 1e-3
    g = lambda f: d["eft_" + f]
    D = lambda x: np.gradient(x, la, edge_order=2)  # d/dln a
    fd = {
        "Hdotdotdot": H * D(g("Hdotdot")),
        "EFTcdotdot": H * D(g("EFTcdot")) - 2 * H * g("EFTcdot"),
        "EFTLambdadotdot": H * D(g("EFTLambdadot")) - 2 * H * g("EFTLambdadot"),
        "EFTGamma1PP": D(g("EFTGamma1P")) / a,
        "EFTOmegaPPPP": D(g("EFTOmegaPPP")) / a,
    }
    out = []
    for f, n in fd.items():
        c = g(f)[s]; n = n[s]
        mc, mn = np.abs(c).max(), np.abs(n).max()
        if mn == 0 and mc == 0:
            out.append(f"{f}=both0"); continue
        rel = np.abs(c - n).max() / max(mn, 1e-300)
        tag = "ZERO" if mc == 0 else ""
        out.append(f"{f}:{tag}rel={rel:.2g}")
    print(f"{m:26s} " + "  ".join(out))

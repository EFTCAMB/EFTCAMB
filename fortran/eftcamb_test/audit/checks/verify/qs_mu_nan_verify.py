"""Independent check of qs-mu-nan-vanishing-denominator: harness get_scale_evolution path (evol_k) on a dense a-grid,
plus a sympy proof that the closed-form QS_mu / QS_mu_gen are 0/0 when OmegaP=c=Gamma2=Gamma5(+derivs)=0."""
import sys, json
import numpy as np
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
A0 = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0"
OUT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/qsnan"
PAR = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/"
models = sys.argv[1:] or ["2_PEFT_gamma1_const_1", "2_PEFT_gamma6_const_1", "5_ADE_1", "1_GR"]
for m in models:
    ini = PAR + m + ".ini"
    npz = f"{OUT}/{m}.npz"
    st = harness.run_model(A0, ini, npz, evol_k=[1e-3, 0.01, 0.1, 1.0], threads=1, timeout=900)
    if st["status"] != "ok":
        print(m, st["status"], st.get("log_tail", "")[-500:]); continue
    d = harness.load(npz)
    ea = d["evol_a"]
    act = ea > 1e-6
    print(m, "EFTflag", (d["meta"].get("EFTflag") if isinstance(d.get("meta"), dict) else "?"))
    for f in ["QS_mu", "QS_sigma", "QS_mu_gen", "QS_sigma_gen", "mu", "gamma", "sigma_eft", "EFTOmegaP", "EFTc"]:
        key = "evol_" + f
        if key not in d: print("  missing", key); continue
        v = np.asarray(d[key])
        nf = ~np.isfinite(v)
        fin = v[np.isfinite(v)]
        print(f"  {f:13s} nonfinite {nf.sum()}/{v.size}  finite-range [{fin.min() if fin.size else np.nan:.4g}, {fin.max() if fin.size else np.nan:.4g}]"
              + (f"  first bad (k_idx,a)={[(float(d["evol_k"][i, j]), float(ea[i, j])) for i, j in zip(*np.where(nf))][:3]}" if nf.any() else ""))

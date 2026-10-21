"""Independent check of scg-adotoa-drops-massive-nu: rebuild H from CAMB's own background densities
(get_background_densities, includes massive nu) + EFT functions, compare to EFTCAMB adotoa; and check that the
deficit in grhom_t equals exactly the massive-nu density."""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(1, os.path.join(os.path.dirname(__file__), ".."))
import harness
import numpy as np
TREE, INI = sys.argv[1], sys.argv[2]
with harness.slot(1):
    sys.path.insert(0, TREE)
    import camb
    cp = camb.read_ini(INI, no_validate=True)
    if not hasattr(cp, "EFTCAMB") or cp.EFTCAMB is None:
        print("no EFTCAMB"); sys.exit()
    res = camb.get_background(cp)
    a = np.array([1e-7, 1e-5, 1e-3, 1e-2, 0.1, 0.5, 1.0])
    _, arr = cp.EFTCAMB.get_eft_functions(res, a)
    d = res.get_background_densities(a, vars=["tot", "nu", "de", "photon", "neutrino", "baryon", "cdm"], format="dict")
# d[x] = 8 pi G rho a^4 (Mpc^-2)
nu_a2 = d["nu"] / a**2  # massive nu, in grhom_t units (8piG rho a^2)
mat_a2 = (d["photon"] + d["neutrino"] + d["baryon"] + d["cdm"] + d["nu"]) / a**2
eft_noNu = mat_a2 - nu_a2
print("a  grhom_t/(CAMB matter w/o massive nu)-1   grhom_t/(CAMB matter incl nu)-1   adotoa_EFT/H_indep-1")
for i in range(len(a)):
    OV, OP, c, L = arr["EFTOmegaV"][i], arr["EFTOmegaP"][i], arr["EFTc"][i], arr["EFTLambda"][i]
    Hind = np.sqrt((mat_a2[i] + 2 * c - L) / 3 / (1 + OV + a[i] * OP))
    print(f"{a[i]:8.1e}  {arr['grhom_t'][i]/eft_noNu[i]-1:+.2e}   {arr['grhom_t'][i]/mat_a2[i]-1:+.2e}   "
          f"{arr['adotoa'][i]/Hind-1:+.3e}")

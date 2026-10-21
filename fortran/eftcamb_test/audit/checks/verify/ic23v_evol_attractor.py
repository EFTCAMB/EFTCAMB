"""Physics-lens check for the ic23 finding: along the evolved solution (get_scale_evolution), compare pi with the
quasi-static attractor -h0 E/(C + k^2 D1 + k^4 D2) and with the case-2 radiation series h0/4 (1+O)/(5+3O) k^2 tau^3.
usage: python ic23v_evol_attractor.py run.npz [run2.npz ...]"""

import sys

import numpy as np

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

for path in sys.argv[1:]:
    d = harness.load(path)
    ks, a = d["evol_k"][:, 0], d["evol_a"][0]
    H0 = 70.0  # base_params.ini hubble = 70
    h0 = H0 / 299792.458
    print(f"== {path}  H0={H0}")
    keys = [k for k in d if k.startswith("evol_")]
    need = ["evol_pi", "evol_EFTpiC", "evol_EFTpiD1", "evol_EFTpiD2", "evol_EFTpiE", "evol_tau", "evol_EFTOmegaV"]
    miss = [n for n in need if n not in d]
    if miss:
        print("missing", miss, "have", sorted(keys)[:40])
        continue
    for i, k in enumerate(ks):
        pi = d["evol_pi"][i]
        C, D1, D2, E = d["evol_EFTpiC"][i], d["evol_EFTpiD1"][i], d["evol_EFTpiD2"][i], d["evol_EFTpiE"][i]
        tau = d["evol_tau"][i]
        Om = d["evol_EFTOmegaV"][i]
        with np.errstate(divide="ignore", invalid="ignore"):
            att = -h0 * E / (C + k**2 * D1 + k**4 * D2)
            ser = h0 * 0.25 * (1 + Om) / (5 + 3 * Om) * k**2 * tau**3
        sel = [np.argmin(abs(a - x)) for x in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0)]
        print(f" k={k:g}")
        for j in sel:
            print(f"   a={a[j]:.2e} tau={tau[j]:9.4g} k*tau={k*tau[j]:9.3g} pi={pi[j]: .4e} attr={att[j]: .4e} "
                  f"series={ser[j]: .4e} pi/attr={pi[j]/att[j] if att[j] else np.nan: .3g}")

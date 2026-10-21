"""Independent check of finding eft-background-densities-bypass.

Claim: CAMBdata.get_background_densities / get_Omega / get_dark_energy_rho_w ignore the EFT background when EFTflag != 0.

Independent methods (different code paths from the repro, which compares to hubble_parameter):
  (1) H from distances: H = c / (d chi / dz), chi from comoving_radial_distance (Romberg integral of dtauda).
  (2) EFT effective fluid from EFTCAMB's own getter get_eft_functions: 8 pi G rho_Q a^4 = grhoq a^2, w_Q = gpresq/grhoq,
      plus Friedmann closure 3 adotoa^2 = grhom_t + grhoq inside EFTCAMB.
  (3) get_Omega('cdm') vs the exact omch2/h^2 (1+z)^3 (H0/H)^2.
Runs on A0 (code under audit) and B1 (pre-rebase EFTCAMB) to attribute the origin.

Usage: python verify_bgdens_bypass.py [A0|B1] [ini ...]
"""

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import harness  # noqa: E402

PARAMS = str(harness.PARAMS_DIR)

WORKER = r"""
import sys, json
import numpy as np
tree, params, checks = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, tree)
sys.path.insert(1, checks)
import camb
from harness import parse_ini
c = 299792.458
z = np.array([0.0, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 100.0])
a = 1 / (1 + z)
out = {"camb": camb.__file__}
for ini in sys.argv[4:]:
    path = params + "/" + ini
    cp = camb.read_ini(path, no_validate=True)
    keys = parse_ini(path)
    fallback = False
    if cp.EFTCAMB.EFTflag == 0 and int(keys.get("EFTflag", "0")) != 0:
        cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)
        fallback = True
    res = camb.get_background(cp)
    H_api = np.asarray(res.hubble_parameter(z))
    # (1) H from distances: central difference of chi(z)
    dz = 1e-4 * (1 + z)
    chi_p = np.asarray(res.comoving_radial_distance(z + dz))
    chi_m = np.asarray(res.comoving_radial_distance(np.maximum(z - dz, 0)))
    H_dist = c * (z + dz - np.maximum(z - dz, 0)) / (chi_p - chi_m)
    d = res.get_background_densities(a)
    H_dens = np.sqrt(d["tot"] / 3) / a**2 * c
    rho_api, w_api = res.get_dark_energy_rho_w(a)
    om_cdm_api = np.asarray(res.get_Omega("cdm", z))
    h = cp.H0 / 100
    om_cdm_true = cp.omch2 / h**2 * (1 + z) ** 3 * (cp.H0 / H_api) ** 2
    om_de_api = np.asarray(res.get_Omega("de", z))
    r = dict(EFTflag=int(cp.EFTCAMB.EFTflag), fallback=fallback, z=z.tolist(), H_api=H_api.tolist(),
             H_dist=H_dist.tolist(), H_dens=H_dens.tolist(), de_api=d["de"].tolist(), w_api=np.asarray(w_api).tolist(),
             rho_api=np.asarray(rho_api).tolist(), om_cdm_api=om_cdm_api.tolist(), om_cdm_true=om_cdm_true.tolist(),
             om_de_api=om_de_api.tolist(), tot_minus_de=(d["tot"] - d["de"]).tolist(),
             de_type=type(cp.DarkEnergy).__name__, de_w=float(getattr(cp.DarkEnergy, "w", float("nan"))))
    if cp.EFTCAMB.EFTflag != 0:
        fields, arr = cp.EFTCAMB.get_eft_functions(res, a)
        f = {k: np.asarray(arr[k]) for k in ("grhoq", "gpresq", "adotoa", "grhom_t")}
        r["de_eft"] = (f["grhoq"] * a**2).tolist()  # 8 pi G rho_Q a^4
        r["w_eft"] = (f["gpresq"] / f["grhoq"]).tolist()
        r["H_eft"] = (f["adotoa"] / a * c).tolist()
        r["friedmann_closure"] = (3 * f["adotoa"] ** 2 / (f["grhom_t"] + f["grhoq"]) - 1).tolist()
        # DE implied by dtauda: 3 (H a^2/c)^2 - (non-DE sum from the same API call)
        r["de_from_H"] = (3 * (H_api / c) ** 2 * a**4 - (d["tot"] - d["de"])).tolist()
    out[ini] = r
print("@@" + json.dumps(out))
"""


def main():
    tree_name = sys.argv[1] if len(sys.argv) > 1 else "A0"
    inis = sys.argv[2:] or ["1_EFT_GR.ini", "2_PEFT_wDE_wCDM_1.ini", "2_PEFT_Omega_const_1.ini",
                            "4_Des_F_R_LCDM_1.ini", "4_Des_mc5e_CPL_1.ini", "3_RPH_KT_OmegaDE_1.ini",
                            "5_Kmouflage_1.ini", "5_quint_exponential.ini", "5_ADE_1.ini"]
    tree = str(harness.TREES[tree_name])
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-c", WORKER, tree, PARAMS, os.path.join(HERE, "..")] + inis,
                           capture_output=True, text=True, env=env, preexec_fn=lambda: os.nice(10), timeout=1500)
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@")]
    if not line:
        print(p.stdout[-3000:], p.stderr[-3000:])
        sys.exit(1)
    res = json.loads(line[0][2:])
    print("tree", tree_name, res.pop("camb"))
    for ini, r in res.items():
        import numpy as np

        H_api, H_dist, H_dens = (np.array(r[k]) for k in ("H_api", "H_dist", "H_dens"))
        print(f"\n{ini}  EFTflag={r['EFTflag']} fallback={r['fallback']} DarkEnergy={r['de_type']} w={r['de_w']}")
        print(f"   max|H_dist/H_api-1| = {np.max(abs(H_dist / H_api - 1)):.2e}   (distance path vs hubble_parameter)")
        print(f"   max|H_dens/H_dist-1| = {np.max(abs(H_dens / H_dist - 1)):.2e}  (get_background_densities vs distances)")
        print(f"   max|Om_cdm_api/Om_cdm_true-1| = {np.max(abs(np.array(r['om_cdm_api']) / np.array(r['om_cdm_true']) - 1)):.2e}")
        print(f"   w_api range = [{min(r['w_api']):.4f}, {max(r['w_api']):.4f}]")
        if "de_eft" in r:
            de_api, de_eft, de_H = (np.array(r[k]) for k in ("de_api", "de_eft", "de_from_H"))
            print(f"   max|H_eft/H_api-1| = {np.max(abs(np.array(r['H_eft']) / H_api - 1)):.2e}   (EFTCAMB adotoa vs dtauda)")
            print(f"   max|friedmann closure| = {np.max(abs(np.array(r['friedmann_closure']))):.2e}")
            low = np.array(r["z"]) <= 10
            print("   z           ", r["z"])
            print("   de_api/de_eft  ", ["%.4f" % x for x in de_api / de_eft])
            print("   de_H/de_eft    ", ["%.4f" % x for x in de_H / de_eft])
            print("   w_api          ", ["%.4f" % x for x in r["w_api"]])
            print("   w_eft          ", ["%.4f" % x for x in r["w_eft"]])
            print(f"   max_z<=10 |de_api/de_eft-1| = {np.max(abs(de_api[low] / de_eft[low] - 1)):.3e}; "
                  f"max |w_api-w_eft| = {np.max(abs(np.array(r['w_api']) - np.array(r['w_eft']))):.3e}")
            print("   Om_de_api      ", ["%.4f" % x for x in r["om_de_api"]])
            print("   Om_de_true     ", ["%.4f" % x for x in de_H / (3 * (H_api / 299792.458) ** 2 * (1 / (1 + np.array(r['z']))) ** 4)])


if __name__ == "__main__":
    main()

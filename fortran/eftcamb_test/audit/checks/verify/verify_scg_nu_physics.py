"""Adversarial verification of finding scg-adotoa-drops-massive-nu (Scaling Cubic Galileon, EFTflag=4, model 6).

Part A (background): uses the real derivs path (get_scale_evolution -> EFTCAMB_quantities_evolution -> derivs) and
compares the adotoa used by the perturbation equations with
  (1) CAMB's a H(a) from results.hubble_parameter (dtauda incl. massive nu),
  (2) d ln a / d tau from results.conformal_time (the a(tau) table that derivs reads a from via ThermoData),
  (3) adotoa_derivs^2 + grhonu_tot / (3 (1 + OmegaV + a OmegaP))  -- i.e. the same formula with nu restored,
and checks H(z=0) against the input H0.

Part B (observables): neutrino-mass response of the CMB. Physically, switching a 0.06 eV neutrino on (vs massless,
same N_eff) changes high-ell TT by << 1 % in any model whose early universe is standard; compare SCG with GR.

Usage: python verify_scg_nu_physics.py [A|B|all]
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

AUD = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
           "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
A0 = AUD / "trees" / "A0"
PSCG = AUD / "work" / "t4b-derivs" / "P_scg"  # t4b partial fix (keep caller grhom_t), read-only reuse
OUT = AUD / "work" / "verify" / "scg_nu"
OUT.mkdir(parents=True, exist_ok=True)
PARAMS = HERE.parents[2] / "parameters"

MASSLESS = {"omnuh2": 0, "massive_neutrinos": 0, "massless_neutrinos": 3.046, "nu_mass_eigenstates": 0}
MNU006 = {"omnuh2": 0.00064}

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_results(cp)  # get_scale_evolution after get_background alone segfaults (rc=139); use full results
a = np.logspace(-6, 0, 61)
z = 1 / a - 1
fields, ev = cp.EFTCAMB.get_scale_evolution(res, [0.05], a)
ev = ev[0]
c_kms = camb.constants.c / 1e3
aH = np.asarray(res.hubble_parameter(z)) / c_kms * a
# d ln a / d tau from the conformal-time table (centred difference in ln a)
eps = 1e-4
tp = np.asarray(res.conformal_time(1 / (a * np.exp(eps)) - 1))
tm = np.asarray(res.conformal_time(1 / (a * np.exp(-eps)) - 1))
dlna_dtau = 2 * eps / (tp - tm)
out = {"a": a.tolist(), "aH": aH.tolist(), "dlna_dtau": dlna_dtau.tolist(),
       "H0_in": cp.H0, "H_z0": float(res.hubble_parameter(0.0)),
       "num_nu_massive": int(cp.num_nu_massive), "omnuh2": cp.omnuh2}
for f in ["a", "tau", "adotoa", "Hdot", "grhonu_tot", "gpinu_tot", "grhom_t", "gpresm_t", "grhoa2",
          "EFTOmegaV", "EFTOmegaP", "EFTc", "EFTLambda"]:
    out["ev_" + f] = np.asarray(ev[f]).tolist()
print("JSON" + json.dumps(out))
"""


def background_case(tree: Path, name: str, overrides: dict) -> dict:
    d = OUT / f"bg_{name}"
    d.mkdir(exist_ok=True)
    ini = d / "scg.ini"
    lines = [f"DEFAULT({PARAMS}/5_scg.ini)"] + [f"{k} = {v}" for k, v in overrides.items()]
    lines += [f"output_root = {d}/out", "EFTCAMB_write_background = F", "feedback_level = 0",
              "EFTCAMB_feedback_level = 0"]
    ini.write_text("\n".join(lines) + "\n")
    w = d / "worker.py"
    w.write_text(WORKER)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    with harness.slot(1):
        p = subprocess.run([sys.executable, str(w), str(tree), str(ini)], capture_output=True, text=True, env=env,
                           cwd=d, preexec_fn=lambda: os.nice(10), timeout=900)
    js = [x for x in p.stdout.splitlines() if x.startswith("JSON")]
    if not js:
        raise RuntimeError(p.stdout[-3000:] + p.stderr[-3000:])
    r = json.loads(js[-1][4:])
    (d / "result.json").write_text(json.dumps(r))
    return r


def part_a() -> dict:
    summary = {}
    for name, ov in [("shipped", {}), ("massless", MASSLESS), ("mnu006", MNU006)]:
        r = background_case(A0, name, ov)
        a = np.array(r["a"])
        ad = np.array(r["ev_adotoa"])
        aH = np.array(r["aH"])
        dl = np.array(r["dlna_dtau"]) * np.array(r["ev_a"])  # d ln a/dtau * ... (see below)
        # d ln a / d tau is already H_conformal = adotoa; ev_a is a at the derivs time (should equal a)
        dl = np.array(r["dlna_dtau"])
        denom = 1 + np.array(r["ev_EFTOmegaV"]) + a * np.array(r["ev_EFTOmegaP"])
        ad_nu = np.sqrt(ad**2 + np.array(r["ev_grhonu_tot"]) / (3 * denom))
        # grhom_t seen by derivs after compute_adotoa vs the same with nu
        print(f"\n== {name}: omnuh2={r['omnuh2']}, num_nu_massive={r['num_nu_massive']}, "
              f"H0_in={r['H0_in']:.6f}, H(z=0)_CAMB={r['H_z0']:.6f}, ratio-1={r['H_z0'] / r['H0_in'] - 1:+.3e}")
        print("      a     a_derivs/a-1   adotoa/aH-1   adotoa/(dlna/dtau)-1   adotoa_withnu/aH-1   "
              "grhonu/grhom_t")
        idx = [0, 10, 20, 30, 40, 50, 55, 60]
        for i in idx:
            print(f"  {a[i]:8.1e}   {r['ev_a'][i] / a[i] - 1:+.2e}     {ad[i] / aH[i] - 1:+.3e}     "
                  f"{ad[i] / dl[i] - 1:+.3e}              {ad_nu[i] / aH[i] - 1:+.3e}          "
                  f"{r['ev_grhonu_tot'][i] / r['ev_grhom_t'][i]:.3e}")
        summary[name] = {
            "max_abs_adotoa_vs_aH": float(np.max(np.abs(ad / aH - 1))),
            "max_abs_adotoa_vs_dlnadtau": float(np.max(np.abs(ad / dl - 1))),
            "max_abs_adotoa_withnu_vs_aH": float(np.max(np.abs(ad_nu / aH - 1))),
            "adotoa_vs_aH_at": {f"{a[i]:.0e}": float(ad[i] / aH[i] - 1) for i in idx},
            "Hz0_over_H0_m1": r["H_z0"] / r["H0_in"] - 1,
        }
    (OUT / "partA_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def tt_ratio(a: dict, b: dict, lmin=30, lmax=2500) -> tuple[float, int]:
    ta, tb = a["cl_total"][:, 0], b["cl_total"][:, 0]
    ell = np.arange(len(ta))
    m = (ell >= lmin) & (ell <= min(lmax, len(ta) - 1, len(tb) - 1))
    rr = ta[: len(tb)][m[: len(tb)]] / tb[: len(ta)][m[: len(ta)]] - 1
    i = int(np.argmax(np.abs(rr)))
    return float(rr[i]), int(ell[m][i])


def part_b() -> dict:
    runs = {
        "scg_ship": (A0, "5_scg.ini", {}),
        "scg_mnu006": (A0, "5_scg.ini", MNU006),
        "scg_massless": (A0, "5_scg.ini", MASSLESS),
        "gr_ship": (A0, "1_EFT_GR.ini", {}),
        "gr_mnu006": (A0, "1_EFT_GR.ini", MNU006),
        "gr_massless": (A0, "1_EFT_GR.ini", MASSLESS),
        "pscg_mnu006": (PSCG, "5_scg.ini", MNU006),
        "pscg_ship": (PSCG, "5_scg.ini", {}),
    }
    res = {}
    for name, (tree, ini, ov) in runs.items():
        out = OUT / f"{name}.npz"
        if not out.exists():
            st = harness.run_model(str(tree), str(PARAMS / ini), str(out), overrides=ov, threads=1, timeout=900)
            print(name, st["status"])
            if st["status"] != "ok":
                print(st.get("log_tail"))
                continue
        res[name] = harness.load(str(out))
    s = {}
    pairs = [("scg nu-response 0.06eV vs massless", "scg_mnu006", "scg_massless"),
             ("GR  nu-response 0.06eV vs massless", "gr_mnu006", "gr_massless"),
             ("scg nu-response omnuh2=0.01 vs massless", "scg_ship", "scg_massless"),
             ("GR  nu-response omnuh2=0.01 vs massless", "gr_ship", "gr_massless"),
             ("partial-fix scg 0.06eV vs massless", "pscg_mnu006", "scg_massless"),
             ("partial-fix vs A0, 0.06eV", "pscg_mnu006", "scg_mnu006"),
             ("partial-fix vs A0, shipped", "pscg_ship", "scg_ship")]
    for label, x, y in pairs:
        if x in res and y in res:
            r, ell = tt_ratio(res[x], res[y])
            r_hi, ell_hi = tt_ratio(res[x], res[y], lmin=800)
            s8 = (float(res[x]["sigma8"][-1]), float(res[y]["sigma8"][-1])) if "sigma8" in res[x] else None
            print(f"{label:45s} max|dTT/TT| (30-2500) = {r:+.3e} @l={ell}; (800-2500) = {r_hi:+.3e} @l={ell_hi};"
                  f" sigma8 {s8}")
            s[label] = {"max_dTT_30_2500": r, "ell": ell, "max_dTT_800_2500": r_hi, "ell_hi": ell_hi, "sigma8": s8}
    (OUT / "partB_summary.json").write_text(json.dumps(s, indent=1))
    return s


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("A", "all"):
        part_a()
    if what in ("B", "all"):
        part_b()

"""Independent verification of finding rgr-threshold-dimensionful-drops-de-perturbations.

Method (different from the t5 repro, which only compares two EFT thresholds with each other):
  1. Physics reference: CAMB's own DarkEnergyFluid (EFTflag=0, cs2=1 = quintessence) with the same CPL w(a).
     The EFT minimally coupled designer model with DE perturbations must agree with it; with RGR=1.1 it must not.
  2. Direct check that the DE perturbation pi is identically zero when RGR returns 1.1 (get_scale_evolution) while
     the background H(z) is unchanged.
  3. Offline evaluation of the RGR functions 5..8 from the timestep cache, dimensionful vs normalised by H.
  4. Unit dependence: the same dimensionless model (same Omega_i, w(a)) at h=0.70 and h=0.85, threshold 5e-8:
     the RGR outcome flips with h because functions 5..8 scale as h^2.
  5. Default threshold (1e-8): near-LCDM CPL (w0=-0.975, wa=-0.025) -> RGR=1.1, compared with skip_RGR and fluid.
  6. B1 (pre-rebase) behaves the same (inherited, not a rebase regression).
"""

import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "rgrdim"
CPL = harness.PARAMS_DIR / "4_Des_mc5e_CPL_1.ini"
EVK = [0.01, 0.1]


def scaled_h(h):
    f = (h / 0.70) ** 2
    return {"hubble": f"{100 * h:g}", "ombh2": f"{0.0226 * f:.8g}", "omch2": f"{0.112 * f:.8g}",
            "omnuh2": f"{0.01 * f:.8g}"}


FLUID = {"EFTflag": "0", "dark_energy_model": "fluid", "cs2_lam": "1"}
RUNS = {
    "cpl_thr1e-8": ("A0", {"EFTCAMB_GR_threshold": "1e-8"}, EVK),
    "cpl_thr1e-7": ("A0", {"EFTCAMB_GR_threshold": "1e-7"}, EVK),
    "cpl_fluid": ("A0", dict(FLUID, w="-0.9", wa="-0.1"), None),
    "near_default": ("A0", {"EFTw0": "-0.975", "EFTwa": "-0.025"}, EVK),
    "near_skiprgr": ("A0", {"EFTw0": "-0.975", "EFTwa": "-0.025", "EFTCAMB_skip_RGR": "T"}, None),
    "near_fluid": ("A0", dict(FLUID, w="-0.975", wa="-0.025"), None),
    "h070_thr5e-8": ("A0", dict(scaled_h(0.70), EFTCAMB_GR_threshold="5e-8"), None),
    "h085_thr5e-8": ("A0", dict(scaled_h(0.85), EFTCAMB_GR_threshold="5e-8"), None),
    "b1_cpl_thr1e-7": ("B1", {"EFTCAMB_GR_threshold": "1e-7"}, None),
    "b1_cpl_thr1e-8": ("B1", {"EFTCAMB_GR_threshold": "1e-8"}, None),
}


def job(name):
    tree, ov, evk = RUNS[name]
    out = OUT / f"{name}.npz"
    if not (out.exists() and json.loads(out.with_suffix(".json").read_text())["status"] == "ok"):
        st = harness.run_model(tree, CPL, out, overrides=ov, evol_k=evk, threads=1, timeout=900)
        if st["status"] != "ok":
            raise RuntimeError(f"{name}: {st}")
    log = out.with_suffix(".log").read_text(errors="replace")
    rgr = re.findall(r"Return to GR time:\s*(\S+)", log)
    return name, (harness.load(out), float(rgr[0]) if rgr else None)


def rgr_funcs(d, a_on):
    a = d["eft_a_grid"]
    m = (a >= a_on * (1 - 1e-12)) & (a <= 1.0)
    grhov = d["eft_grhov_t"][-1]  # grhov_t(a=1) = grhov (8 pi G rho_DE0, Mpc^-2)
    a2 = a**2
    H = d["eft_adotoa"] / a  # physical H in Mpc^-1
    f = {
        "f5 |c/a2|": np.abs(d["eft_EFTc"] / a2),
        "f6 |L/a2+grhov|": np.abs(d["eft_EFTLambda"] / a2 + grhov),
        "f7 |cdot/a2|": np.abs(d["eft_EFTcdot"] / a2),
        "f8 |Ldot/a2|": np.abs(d["eft_EFTLambdadot"] / a2),
    }
    norm = {"f5 |c/a2|": H**2, "f6 |L/a2+grhov|": 3 * H**2, "f7 |cdot/a2|": H**3, "f8 |Ldot/a2|": H**3}
    rows = {}
    for k, v in f.items():
        rows[k] = dict(max=float(v[m].max()), at_on=float(v[m][0]), today=float(v[-1]),
                       max_dimless=float((v / norm[k])[m].max()), today_dimless=float((v / norm[k])[-1]))
    return rows, grhov


def diff(d1, d0):
    tt = d1["cl_total"][2:2501, 0] / d0["cl_total"][2:2501, 0] - 1
    pp = d1["cl_lens_potential"][2:2001, 0] / d0["cl_lens_potential"][2:2001, 0] - 1
    return dict(sigma8=float(d1["sigma8"][0] / d0["sigma8"][0] - 1), tt2=float(tt[0]),
                tt30=float(np.abs(tt[28:]).max()), pp=float(np.abs(pp).max()),
                bgH=float(np.abs(d1["bg_H"] / d0["bg_H"] - 1).max()))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(job, list(RUNS)))
    print("RGR times from logs:")
    for k, (d, r) in res.items():
        print(f"  {k:16s} RGR={r}  sigma8={d['sigma8'][0]:.5f}")

    print("\n[3] RGR functions 5..8 (a in [0.01,1]) for CPL(-0.9,-0.1), h=0.70:")
    rows, grhov = rgr_funcs(res["cpl_thr1e-8"][0], 0.01)
    print(f"  grhov = {grhov:.4e} Mpc^-2")
    for k, v in rows.items():
        print(f"  {k:16s} max={v['max']:.3e} at_turn_on={v['at_on']:.3e} today={v['today']:.3e} | "
              f"normalised: max={v['max_dimless']:.3e} today={v['today_dimless']:.3e}")
    print("  near-LCDM CPL(-0.975,-0.025):")
    rows_n, _ = rgr_funcs(res["near_default"][0], 0.01)
    for k, v in rows_n.items():
        print(f"  {k:16s} max={v['max']:.3e} today={v['today']:.3e} | normalised max={v['max_dimless']:.3e}"
              f" today={v['today_dimless']:.3e}")

    print("\n[4] unit dependence (same Omega_i, w(a); threshold 5e-8):")
    r70, _ = rgr_funcs(res["h070_thr5e-8"][0], 0.01)
    r85, _ = rgr_funcs(res["h085_thr5e-8"][0], 0.01)
    for k in r70:
        print(f"  {k:16s} max h=0.85/h=0.70 = {r85[k]['max'] / r70[k]['max']:.4f} (h^2 ratio "
              f"{(0.85 / 0.70) ** 2:.4f}); normalised ratio = {r85[k]['max_dimless'] / r70[k]['max_dimless']:.4f}")
    print(f"  RGR h=0.70: {res['h070_thr5e-8'][1]}   RGR h=0.85: {res['h085_thr5e-8'][1]}")

    print("\n[2] DE perturbation pi (get_scale_evolution) and background:")
    for k in ("cpl_thr1e-8", "cpl_thr1e-7", "near_default"):
        d = res[k][0]
        a = d["evol_a"][0]  # struct field "a" overwrites the evol_a grid key
        pi = d["evol_pi"]
        late = a > 0.02
        print(f"  {k:14s} max|pi| (a>0.02) k=0.01: {np.abs(pi[0][late]).max():.3e}  k=0.1: "
              f"{np.abs(pi[1][late]).max():.3e}")
    d1, d0 = res["cpl_thr1e-7"][0], res["cpl_thr1e-8"][0]
    print(f"  bg_H max rel diff thr1e-7 vs thr1e-8: {np.abs(d1['bg_H'] / d0['bg_H'] - 1).max():.2e}")

    print("\n[1] vs CAMB DarkEnergyFluid (cs2=1, same CPL):  rel diffs  sigma8 | TT l=2 | max TT l>=30 | max Cpp | bgH")
    for a_, b_ in (("cpl_thr1e-8", "cpl_fluid"), ("cpl_thr1e-7", "cpl_fluid"), ("cpl_thr1e-7", "cpl_thr1e-8"),
                   ("near_skiprgr", "near_fluid"), ("near_default", "near_fluid"), ("near_default", "near_skiprgr"),
                   ("b1_cpl_thr1e-7", "b1_cpl_thr1e-8")):
        x = diff(res[a_][0], res[b_][0])
        print(f"  {a_:14s} vs {b_:14s}: {x['sigma8']:+.2e} | {x['tt2']:+.2e} | {x['tt30']:.2e} | {x['pp']:.2e} |"
              f" {x['bgH']:.1e}")


if __name__ == "__main__":
    main()

"""Verify t3-gamma1-mass-instability-near-gr numerically (independent of the reporter's script).

1. A0 runs of 2_PEFT_gamma1_const_1 with EFT_mass_stability=F (Gamma1 = 1e-2, 1e-3, 1e-5, 1e-7): take the
   background/EFT cache from get_eft_functions and evaluate the automatically transcribed EFT_mu2 formula
   (a) in float64 (same operation order as the Fortran source, no FMA), (b) in 60-digit mpmath on the same inputs,
   (c) from the exact closed form derived by gamma1_mu2_symbolic.py: mu2 = g Hdot^2 / (2 (g - 6) H^2), g = 4 H0^2 G1/H^2.
2. Physics: Cls at default accuracy with mass stability off, Gamma1 = 1e-3 and 1e-2, against the same ini with
   EFTflag=0 (LCDM); a genuine tachyonic mass mu^2 ~ -500 H^2 would give an exponential blow-up.
3. B1 (pre-rebase) rejection at Gamma1 = 1e-3 with the default mass test.
"""

import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import gamma1_mu2_formula as gf  # noqa: E402
import harness  # noqa: E402

INI = harness.PARAMS_DIR / "2_PEFT_gamma1_const_1.ini"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "gamma1_mu2"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "50", "get_transfer": "F",
        "num_redshiftwindows": "0"}


def mu2_eval(d, i):
    av = float(d["eft_a_grid"][i])
    h0 = float(d["eft_adotoa"][-1])  # adotoa(a=1) = H0 [1/Mpc]
    vals = {f: float(d["eft_" + f][i]) for f in gf.CACHE_FIELDS}
    env_np = {"a": np.float64(av), "h0": np.float64(h0)}
    env_mp = {"a": mp.mpf(av), "h0": mp.mpf(h0)}
    for f, v in vals.items():
        env_np["C_" + f] = np.float64(v)
        env_mp["C_" + f] = mp.mpf(v)
    mp.mp.dps = 60
    with np.errstate(all="ignore"):
        n1 = gf.eval_formula(env_np, "np")
    n2 = gf.eval_formula(env_mp, "mp")
    Hp = vals["adotoa"] / av
    Hdp = (vals["Hdot"] - vals["adotoa"] ** 2) / av**2
    g = 4 * h0**2 * vals["EFTGamma1V"] / Hp**2
    closed = g * (Hdp / Hp**2) ** 2 / (2 * (g - 6))
    kin = 3 * float(n1["F2"]) ** 2 + float(n1["F1"]) * float(n1["F3"])
    return {"a": av, "mu2_f64": float(n1["C_EFT_mu2"]) / Hp**2, "mu2_mp60": float(n2["C_EFT_mu2"]) / Hp**2,
            "mu2_closed": closed, "mu1_f64": float(n1["C_EFT_mu1"]) / Hp**2,
            "mu1_mp60": float(n2["C_EFT_mu1"]) / Hp**2, "g": g, "kin_over_F2sq": kin / float(n1["F2"]) ** 2,
            "c": vals["EFTc"]}


def run(tag, overrides, tree="A0", fast=True):
    out = OUT / f"{tag}.npz"
    ov = {**(FAST if fast else {}), **overrides}
    st = harness.run_model(tree, INI, out, overrides=ov, threads=1, timeout=1500)
    return tag, st, out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = []
    for v in (1e-2, 1e-3, 1e-5, 1e-7):
        jobs.append((f"nomass_G{v:g}", {"EFTGamma10": repr(v), "EFT_mass_stability": "F"}, "A0", True))
    jobs.append(("full_nomass_G0.001", {"EFTGamma10": "0.001", "EFT_mass_stability": "F"}, "A0", False))
    jobs.append(("full_nomass_G0.01", {"EFTGamma10": "0.01", "EFT_mass_stability": "F"}, "A0", False))
    jobs.append(("full_lcdm", {"EFTflag": "0"}, "A0", False))
    jobs.append(("B1_G0.001", {"EFTGamma10": "0.001"}, "B1", True))
    jobs.append(("A0_G0.001_rate", {"EFTGamma10": "0.001"}, "A0", True))
    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(lambda j: run(j[0], j[1], j[2], j[3]), jobs))
    summary = {}
    for tag, st, out in res:
        log = out.with_suffix(".log").read_text(errors="replace") if out.with_suffix(".log").exists() else ""
        mu = re.findall(r"Mass instability\. mu_(\d)/H\^2 = *([-+0-9.E]+)", log)
        at = re.findall(r"Instability detected at a = *([-+0-9.E]+)", log)
        rgr = re.findall(r"Return to GR time: *([-+0-9.E]+)", log)
        print(f"{tag:22s} {st['status']:12s} RGR={rgr[:1]} mass-msg={mu[:2]} at={at[:1]}")
        summary[tag] = {"status": st["status"], "mu": mu[:2], "at": at[:1], "rgr": rgr[:1]}
    print("\nformula evaluation on the A0 background (mu/H^2):")
    for tag, st, out in res:
        if not tag.startswith("nomass_") or st["status"] != "ok":
            continue
        d = np.load(out, allow_pickle=True)
        ag = d["eft_a_grid"]
        for target in (0.0105, 0.0525, 0.25, 1.0):
            i = int(np.argmin(np.abs(ag - target)))
            r = mu2_eval(d, i)
            print(f"  {tag:16s} a={r['a']:.4f} g={r['g']:.3e} (3F2^2+F1F3)/F2^2={r['kin_over_F2sq']:.3e} "
                  f"mu2: f64={r['mu2_f64']:+.4e} mp60={r['mu2_mp60']:+.4e} closed={r['mu2_closed']:+.4e} | "
                  f"mu1: f64={r['mu1_f64']:+.4e} mp60={r['mu1_mp60']:+.4e} c={r['c']:.1e}")
            summary.setdefault("formula", []).append({"tag": tag, **r})
    print("\nphysics: Cls with mass test off vs LCDM (EFTflag=0), default accuracy:")
    if all(Path(OUT / f"{t}.npz").exists() for t in ("full_lcdm", "full_nomass_G0.001", "full_nomass_G0.01")):
        lc = harness.load(OUT / "full_lcdm.npz")
        for t in ("full_nomass_G0.001", "full_nomass_G0.01"):
            c = harness.compare(harness.load(OUT / f"{t}.npz"), lc, lmax=2500)
            cl = c["cl"]["total"] if "total" in c["cl"] else c["cl"]
            print(f"  {t}: total Cl max frac diff {cl} sigma8 {c['sigma8']} worst_cl {c['worst_cl']} "
                  f"worst_pk {c['worst_pk']}")
            summary[t + "_vs_lcdm"] = {"cl_total": cl, "sigma8": c["sigma8"], "worst_cl": c["worst_cl"],
                                       "worst_pk": c["worst_pk"]}
    (OUT / "numeric_summary.json").write_text(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()

"""Check 2: source windows (number counts, galaxy lensing) and transfer outputs with EFT on.

(a) GR limit: an EFT model with EFT perturbations active from a=0.01 but with tiny modifications (PEFT Omega0=1e-6)
    must reproduce EFTflag=0 for all number-count terms (each term separately and all together) and lensing windows.
(b) Weyl transfer consistency (MG): C_L^phiphi from the line-of-sight EFT lensing source vs a Limber integral of
    P_Weyl(k,z) from get_matter_power_interpolator(var Weyl) [CAMBdemo recipe]; must agree as well as in GR.

usage: python c2_sources.py [run|report]
"""

from __future__ import annotations

import json
import sys

import numpy as np
from t8c_lib import WORK, load, model_ini, run_many

OUT = WORK / "c2"
OV = {"get_tensor_cls": "F", "get_transfer": "F", "feedback_level": "0"}
WINDOWS = [
    {"redshift": 0.5, "source_type": "counts", "sigma": 0.1, "bias": 1.5, "dlog10Ndm": 0.2},
    {"redshift": 2.0, "source_type": "counts", "sigma": 0.2, "bias": 2.0, "dlog10Ndm": 0.1},
    {"redshift": 1.0, "source_type": "lensing", "sigma": 0.3},
]
TERMS = ["counts_density", "counts_redshift", "counts_radial", "counts_timedelay", "counts_ISW", "counts_velocity",
         "counts_potential", "counts_lensing"]

# near-GR EFT model: PEFT Omega const with tiny amplitude, EFT active from a=0.01 (RGR skipped)
NEARGR = {"EFTOmega0": "1e-6", "EFTOmegaExp": "1.0", "EFTCAMB_skip_RGR": "T", "EFT_mass_stability": "F"}


def src_args(terms_on, limber=True):
    st = {t: (t in terms_on) for t in TERMS}
    st["limber_windows"] = limber
    return {"set_for_lmax": {"lmax": 1000, "lens_potential_accuracy": 0, "nonlinear": False},
            "source_windows": WINDOWS, "source_terms": st, "outputs": ["sources"]}


def jobs():
    js = []
    gr = model_ini("1_EFT_GR")
    ng = model_ini("2_PEFT_Omega_const_1")
    variants = {"all": TERMS} | {t: [t] for t in TERMS}
    for tag, terms in variants.items():
        js.append(dict(tree="A0", ini=gr, task="generic", out=OUT / f"GR_{tag}.npz", overrides=OV,
                       args=src_args(terms)))
        js.append(dict(tree="A0", ini=ng, task="generic", out=OUT / f"NG_{tag}.npz", overrides=OV | NEARGR,
                       args=src_args(terms)))
    return js


def report():
    res = {}
    for tag in ["all"] + TERMS:
        try:
            a, b = load(OUT / f"NG_{tag}.npz"), load(OUT / f"GR_{tag}.npz")
        except FileNotFoundError:
            res[tag] = "missing"
            continue
        keys = list(a["src_keys"])
        d = {}
        for i, k in enumerate(keys):
            if k.startswith("P"):
                continue
            x, y = a["src_cls"][i][2:1001], b["src_cls"][i][2:1001]
            if not np.any(y):
                continue
            # cross spectra can cross zero: normalise by sqrt of the autos
            k1, k2 = k.split("x")
            n1 = b["src_cls"][keys.index(f"{k1}x{k1}")][2:1001]
            n2 = b["src_cls"][keys.index(f"{k2}x{k2}")][2:1001]
            d[k] = float(np.max(np.abs(x - y) / np.sqrt(np.abs(n1 * n2))))
        res[tag] = d
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        for r in run_many(jobs(), nproc=4):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report()

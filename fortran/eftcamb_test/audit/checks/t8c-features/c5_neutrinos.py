"""Check 5: several massive neutrino eigenstates and the v2 hierarchy/mass-density setter with EFT on.

set_cosmology(...) is called after read_ini with each neutrino setup (v2: mnu is the physical mass sum, omnuh2 from
the thermal density; omnuh2_active alternative; standard_neutrino_neff applied to degeneracies).
(a) GR limit: near-GR EFT model (PEFT Omega0=1e-6, EFT active from a=0.01) vs EFTflag=0 for each setup.
(b) designer f(R), quintessence: runs, finite, EFT parameter cache neutrino masses == CAMB's.
(c) K-mouflage (no massive-nu support in the model): what happens with the Python default mnu=0.06.
"""

from __future__ import annotations

import json
import sys

import numpy as np
from t8c_lib import WORK, cl_diffs, load, maxrel, model_ini, pk_diff, run_many

OUT = WORK / "c5"
OV = {"get_tensor_cls": "F", "feedback_level": "0", "num_redshiftwindows": "0", "transfer_kmax": "2"}
NEARGR = {"EFTOmega0": "1e-6", "EFTCAMB_skip_RGR": "T", "EFT_mass_stability": "F"}
BASE = {"H0": 70.0, "ombh2": 0.0226, "omch2": 0.112, "YHe": 0.24, "TCMB": 2.7255, "nnu": 3.046}
NU = {
    "deg3": {"mnu": 0.15, "num_massive_neutrinos": 3, "neutrino_hierarchy": "degenerate"},
    "normal": {"mnu": 0.10, "neutrino_hierarchy": "normal"},
    "inverted": {"mnu": 0.12, "neutrino_hierarchy": "inverted"},
    "active": {"mnu": None, "omnuh2_active": 0.002, "neutrino_hierarchy": "normal"},
    "default": {},
}
MODELS = {"GR": ("1_EFT_GR", {}), "NG": ("2_PEFT_Omega_const_1", NEARGR), "fR": ("4_Des_F_R_LCDM_1", {}),
          "quint": ("5_quint_exponential", {})}
OUTS = {"outputs": ["cmb", "pk_lin", "nu"], "cl_lmax": 2500}


def jobs():
    js = []
    for nu, kw in NU.items():
        for tag, (m, ov) in MODELS.items():
            js.append(dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{tag}_{nu}.npz",
                           overrides=OV | ov, args=OUTS | {"cosmology": BASE | kw}))
    js.append(dict(tree="A0", ini=model_ini("5_Kmouflage_1"), task="generic", out=OUT / "Kmouflage_default.npz",
                   overrides=OV, args=OUTS | {"cosmology": BASE}))
    return js


def report():
    res = {}
    for nu in NU:
        row = {}
        try:
            g = load(OUT / f"GR_{nu}.npz")
        except FileNotFoundError:
            res[nu] = "missing"
            continue
        row["GR_nu"] = g["meta"].get("nu")
        for tag in ("NG", "fR", "quint"):
            try:
                a = load(OUT / f"{tag}_{nu}.npz")
            except FileNotFoundError:
                st = json.loads((OUT / f"{tag}_{nu}.json").read_text())
                row[tag] = st["status"]
                continue
            nm = a["meta"].get("nu", {})
            ent = {"finite": bool(np.all(np.isfinite(a["cl_total"])) and np.all(np.isfinite(a["pk"]))),
                   "cache_mass_match": nm.get("cache_match"), "sigma8": float(a["sigma8"][0])}
            if tag == "NG":
                ent |= {"cl": cl_diffs(a["cl_total"], g["cl_total"], lmax=2500), "pk": pk_diff(a, g),
                        "sigma8_rel": maxrel(a["sigma8"], g["sigma8"])}
            row[tag] = ent
        res[nu] = row
    st = json.loads((OUT / "Kmouflage_default.json").read_text())
    res["Kmouflage_default_mnu"] = {"status": st["status"], "tail": st.get("log_tail", [])[-3:]}
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        for r in run_many(jobs(), nproc=int(sys.argv[2]) if len(sys.argv) > 2 else 4):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report()

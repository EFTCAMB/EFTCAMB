"""Check 6: omk != 0 with EFT on (supported or refused?).

No EFTCAMB code refuses omk != 0.  GR-limit test: near-GR EFT model (PEFT Omega0=1e-6, EFT perturbations active from
a=0.01, RGR skipped) vs EFTflag=0, at omk = 0, +0.02, -0.02.  If curvature were supported the differences would stay
at the omk=0 level.  Also 2_PEFT_wDE_LCDM_1 (EFT background solver, perturbations GR since RGR time > 1).
"""

from __future__ import annotations

import json
import sys

from t8c_lib import WORK, cl_diffs, load, maxrel, model_ini, pk_diff, run_many

OUT = WORK / "c6"
OV = {"get_tensor_cls": "F", "feedback_level": "0", "num_redshiftwindows": "0", "transfer_kmax": "2"}
NEARGR = {"EFTOmega0": "1e-6", "EFTCAMB_skip_RGR": "T", "EFT_mass_stability": "F"}
OUTS = {"outputs": ["cmb", "pk_lin", "background"], "cl_lmax": 2500}
OMKS = ["0", "0.02", "-0.02"]


def jobs():
    js = []
    for omk in OMKS:
        for tag, m, ov in (("GR", "1_EFT_GR", {}), ("NG", "2_PEFT_Omega_const_1", NEARGR),
                           ("LCDMdes", "2_PEFT_wDE_LCDM_1", {}), ("fR", "4_Des_F_R_LCDM_1", {})):
            js.append(dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{tag}_omk{omk}.npz",
                           overrides=OV | ov | {"omk": omk}, args=OUTS))
    return js


def report():
    res = {}
    for omk in OMKS:
        try:
            g = load(OUT / f"GR_omk{omk}.npz")
        except FileNotFoundError:
            res[omk] = "missing"
            continue
        row = {}
        for tag in ("NG", "LCDMdes"):
            try:
                a = load(OUT / f"{tag}_omk{omk}.npz")
            except FileNotFoundError:
                row[tag] = "missing"
                continue
            row[tag] = {"cl": cl_diffs(a["cl_total"], g["cl_total"], lmax=2500),
                        "pk": pk_diff(a, g), "sigma8": maxrel(a["sigma8"], g["sigma8"]),
                        "H": maxrel(a["bg_H"], g["bg_H"]), "DA": maxrel(a["bg_DA"][1:], g["bg_DA"][1:])}
        try:
            f = load(OUT / f"fR_omk{omk}.npz")
            row["fR_runs"] = True
            row["fR_sigma8"] = float(f["sigma8"][0])
        except FileNotFoundError:
            row["fR_runs"] = False
        res[omk] = row
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        for r in run_many(jobs(), nproc=int(sys.argv[2]) if len(sys.argv) > 2 else 4):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report()

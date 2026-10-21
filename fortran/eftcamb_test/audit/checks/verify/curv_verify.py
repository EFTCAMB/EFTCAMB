"""Independent check of curvature-silently-wrong-with-eft: harness runs (A0, B1), LCDM-limit EFT vs GR at omk=0.02,
comparing CMB TT, background (bg_H) and the EFT cache adotoa vs the GR expectation sqrt((grho+grhok)/3)."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

S = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
OUT = Path(S) / "work/verify/curv"
OUT.mkdir(parents=True, exist_ok=True)
P = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")
res = {}
for tree in ("A0", "B1"):
    for model in ("1_EFT_GR", "2_PEFT_wDE_LCDM_1"):
        for omk in (0.0, 0.02):
            npz = OUT / f"{tree}_{model}_{omk}.npz"
            if not npz.exists():
                st = harness.run_model(f"{S}/trees/{tree}", str(P / f"{model}.ini"), str(npz),
                                       overrides={"omk": omk, "get_tensor_cls": "F", "num_redshiftwindows": 0}, threads=1)
                print(tree, model, omk, st["status"])
            res[(tree, model, omk)] = harness.load(str(npz))
for tree in ("A0", "B1"):
    for omk in (0.0, 0.02):
        g = res[(tree, "1_EFT_GR", omk)]
        e = res[(tree, "2_PEFT_wDE_LCDM_1", omk)]
        tt = e["cl_total"][:, 0] / g["cl_total"][:, 0] - 1
        dH = np.max(np.abs(e["bg_H"] / g["bg_H"] - 1))
        line = f"{tree} omk={omk}: TT l=2 {tt[2]:+.3e} l=10 {tt[10]:+.3e} l=30 {tt[30]:+.3e} l=1000 {tt[1000]:+.3e}"
        line += f" | sigma8 {np.ravel(e['sigma8'])[-1] / np.ravel(g['sigma8'])[-1] - 1:+.2e} | max bg_H rel {dH:.2e}"
        if "eft_adotoa" in e and "eft_adotoa" in g:
            a = e["eft_a_grid"]
            m = a > 1e-3
            line += f" | adotoa(EFT)/adotoa(GR)-1 at a=1: {e['eft_adotoa'][-1] / g['eft_adotoa'][-1] - 1:+.3e}"
            line += f" max(a>1e-3) {np.max(np.abs(e['eft_adotoa'][m] / g['eft_adotoa'][m] - 1)):.3e}"
        print(line)

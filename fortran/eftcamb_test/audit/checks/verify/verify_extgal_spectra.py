"""Verify (adversarial, physics lens) finding extgal-adotoa-stale-grhonu: observable impact.

Runs the Extended Galileon audit inis on A0 (code under audit) and on t2d-static's FIX3 tree (massive-nu density
computed locally in ComputeAdotoa, real bisection tolerance), with default stability checks, for:
  * GCCG1 (s=2, q=0.35) on base_params (omnuh2 = 0.01)
  * GCCG2 (s=1.3, q=0.5) on base_params
  * GCCG1 massless-neutrino control (omnuh2 = 0): FIX3 must agree with A0 to the bisection floor (~4e-7 in H)
and prints harness.compare(A0, FIX3).
"""

from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

SP = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45"
          "/scratchpad/audit")
TREES = {"A0": SP / "trees/A0", "FIX3": SP / "work/t2d-static/FIX3"}
OUT = SP / "work/verify/extgal/spectra"
INIS = REPO / "fortran/eftcamb_test/audit/inis"
CASES = {
    "gccg1_nu": (INIS / "008p7-ext-galileon.ini", {}),
    "gccg2_nu": (INIS / "008p7-ext-galileon_2.ini", {}),
    "gccg1_massless": (INIS / "008p7-ext-galileon.ini",
                       {"omnuh2": "0", "massive_neutrinos": "0", "massless_neutrinos": "3.046"}),
}


def main() -> None:
    jobs = []
    for case, (ini, ov) in CASES.items():
        for t, tree in TREES.items():
            jobs.append((case, t, tree, ini, ov))

    def one(j):
        case, t, tree, ini, ov = j
        return case, t, harness.run_model(tree, ini, OUT / f"{case}_{t}.npz", overrides=ov, threads=1)

    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(one, jobs))
    for case, t, st in res:
        print(case, t, st["status"], st["elapsed"])
    summary = {}
    for case in CASES:
        a, b = OUT / f"{case}_A0.npz", OUT / f"{case}_FIX3.npz"
        if not (a.exists() and b.exists()):
            print(case, "missing output")
            continue
        c = harness.compare(b, a)  # relative to A0
        A, B = harness.load(a), harness.load(b)
        s = {
            "worst_cl": c["worst_cl"], "cl_lensed": c["cl"].get("lensed_scalar"), "cl_unlensed": c["cl"].get(
                "unlensed_scalar"), "pk_z": c.get("pk"), "sigma8": c.get("sigma8"), "bg_H": c.get("bg_H"),
            "bg_tau": c.get("bg_tau"), "eft_adotoa": c.get("eft", {}).get("adotoa"),
            "derived_A0": {k: A["meta"].get("derived", {}).get(k) for k in ("thetastar", "DAstar", "age", "zstar",
                                                                             "rstar")},
            "derived_FIX3": {k: B["meta"].get("derived", {}).get(k) for k in ("thetastar", "DAstar", "age", "zstar",
                                                                               "rstar")},
            "errors_A0": A["meta"].get("errors"), "errors_FIX3": B["meta"].get("errors"),
        }
        summary[case] = s
        print(f"== {case}")
        print(json.dumps(s, indent=1, default=str))
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()

"""Adversarial verifier, finding extgal-adotoa-stale-grhonu: observable impact of the neutrino part of the defect.

Runs the Extended Galileon audit ini on A0 (code under audit) and on t2d-static's FIX3 tree (massive neutrinos computed
locally in ComputeAdotoa + finite bisection tolerance; the a < 0.2 shortcut is unchanged) for omnuh2 = 0.01 (base_params)
and omnuh2 = 0.00064 (one 0.06 eV eigenstate), and compares Cls, P(k), sigma8, H(z), tau(z) with harness.compare.

Usage: python verify_extgal_cl_impact.py OUTDIR
"""

from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REPO = HERE.parents[4]
INI = REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"
SP = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
          "scratchpad/audit")
TREES = {"A0": SP / "trees/A0", "FIX3": SP / "work/t2d-static/FIX3"}
CASES = {"nu0.01": {}, "nu0.00064": {"omnuh2": "0.00064"}}


def main() -> None:
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    jobs = [(t, c) for t in TREES for c in CASES]

    def one(tc):
        t, c = tc
        return harness.run_model(TREES[t], INI, out / f"{t}_{c}.npz", overrides=CASES[c], threads=1, timeout=1500)

    with ThreadPoolExecutor(4) as ex:
        st = list(ex.map(one, jobs))
    for (t, c), s in zip(jobs, st):
        print(t, c, s["status"], s["elapsed"], s.get("log_tail", "")[-3:] if s["status"] != "ok" else "")
    res = {}
    for c in CASES:
        a, b = out / f"A0_{c}.npz", out / f"FIX3_{c}.npz"
        if a.exists() and b.exists():
            cmp_ = harness.compare(a, b)
            res[c] = {k: cmp_[k] for k in ("cl", "worst_cl", "worst_pk", "sigma8", "bg_H", "bg_tau") if k in cmp_}
            res[c]["eft_adotoa"] = cmp_.get("eft", {}).get("adotoa")
            A, B = harness.load(a), harness.load(b)
            res[c]["derived_A0"] = {k: A["meta"].get("derived", {}).get(k) for k in ("thetastar", "rstar", "zstar")}
            res[c]["derived_FIX3"] = {k: B["meta"].get("derived", {}).get(k) for k in ("thetastar", "rstar", "zstar")}
            print(c, json.dumps(res[c], indent=1))
    (out / "cl_impact.json").write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()

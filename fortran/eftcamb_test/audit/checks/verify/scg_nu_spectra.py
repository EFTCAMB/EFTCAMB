"""Differential spectrum test for scg-adotoa-drops-massive-nu (verifier).

Delta_X = spectra(1 massive nu, omnuh2) vs spectra(3.046 massless nu) for X = SCG (5_scg.ini) and X = GR (same ini,
EFTflag=0).  With a correct massive-neutrino treatment Delta_SCG should be of the same size as Delta_GR (a 0.06 eV
neutrino barely changes the primary CMB).  Usage: python scg_nu_spectra.py OUTDIR [omnuh2]
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

PARAMS = HERE.parents[2] / "parameters"


def main() -> None:
    outdir = Path(sys.argv[1])
    omnu = sys.argv[2] if len(sys.argv) > 2 else "0.00064"
    massive = {"omnuh2": omnu}
    massless = {"omnuh2": "0", "massive_neutrinos": "0", "massless_neutrinos": "3.046", "nu_mass_eigenstates": "0"}
    ini = PARAMS / "5_scg.ini"
    jobs = {
        f"scg_massive_{omnu}": dict(massive),
        "scg_massless": dict(massless),
        f"gr_massive_{omnu}": dict(massive, EFTflag="0"),
        "gr_massless": dict(massless, EFTflag="0"),
    }
    from concurrent.futures import ThreadPoolExecutor

    def one(item):
        name, ov = item
        out = outdir / f"{name}.npz"
        if out.exists():
            return name, {"status": "ok", "cached": True}
        return name, harness.run_model("A0", ini, out, overrides=ov, threads=1)

    with ThreadPoolExecutor(4) as ex:
        st = dict(ex.map(one, jobs.items()))
    print(json.dumps({k: v["status"] for k, v in st.items()}))
    d = {k: harness.load(outdir / f"{k}.npz") for k in jobs}
    res = {}
    for lab, m, l0 in (("SCG", f"scg_massive_{omnu}", "scg_massless"), ("GR", f"gr_massive_{omnu}", "gr_massless")):
        a, b = d[m], d[l0]
        cmp = harness.compare(a, b)
        tt_a, tt_b = a["cl_unlensed_scalar"][2:2501, 0], b["cl_unlensed_scalar"][2:2501, 0]
        rel = tt_a / tt_b - 1
        il = int(np.argmax(np.abs(rel)))
        res[lab] = {
            "unlensed_TT_maxrel": float(np.max(np.abs(rel))), "at_l": il + 2,
            "lensed_total": cmp["cl"].get("total", {}),
            "sigma8": [float(a["sigma8"][-1]), float(b["sigma8"][-1])],
            "H0_derived": [a["meta"]["derived"].get("thetastar"), b["meta"]["derived"].get("thetastar")],
        }
        print(f"{lab}: max|TT_unlensed(massive)/TT(massless)-1| = {res[lab]['unlensed_TT_maxrel']:.4e} at l={il + 2};"
              f" lensed {cmp['cl'].get('total')}; sigma8 {res[lab]['sigma8']}; thetastar {res[lab]['H0_derived']}")
    (outdir / f"summary_{omnu}.json").write_text(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()

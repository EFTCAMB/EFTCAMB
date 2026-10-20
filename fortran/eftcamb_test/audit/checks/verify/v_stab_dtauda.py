"""Verify t3-stability-dtauda-without-matter.

1) Rerun the repro scan (OL gamma approaching GR) on A0, A0+skip_stability, R5 (grhoa2 set in stability),
   using an ini copy whose DEFAULT path points at the current repo location.
2) Independent check: from the skip_stability A0 runs, rebuild the full-model dtauda argument from the
   cached EFT functions with and without the matter term, and compare grhoa2 returned by the python bridge
   with grhom_t * a^2.
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

OUT = harness.AUDIT_ROOT / "work" / "verify" / "stab_dtauda"
OUT.mkdir(parents=True, exist_ok=True)
SRC = REPO / "fortran" / "eftcamb_test" / "audit" / "inis" / "007p5-OL-gamma.ini"
INI = OUT / "007p5-OL-gamma.ini"
txt = re.sub(r"DEFAULT\(.*\)", f"DEFAULT({REPO}/fortran/eftcamb_test/parameters/base_params.ini)", SRC.read_text())
INI.write_text(txt)
R5 = harness.AUDIT_ROOT / "work" / "t3-gr" / "R5"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "100", "get_transfer": "F",
        "num_redshiftwindows": "0"}
POINTS = {"Om1e-2": {"OLOmega0": "1e-2"}, "Om1e-4": {"OLOmega0": "1e-4"}, "Om1e-6": {"OLOmega0": "1e-6"},
          "Om0": {"OLOmegamodel": "0"}}


def job(args):
    tree, tag, name, ov = args
    out = OUT / f"{tag}_{name}.npz"
    st = harness.run_model(tree, INI, out, overrides={**FAST, **ov}, threads=1, timeout=900)
    log = out.with_suffix(".log").read_text(errors="replace")
    why = re.findall(r"(Model dtauda is Nan|Physical instability[^\n]*|Mathematical instability[^\n]*|"
                     r"[Gg]host[^\n]*|kinetic[^\n]*)", log)
    at = re.findall(r"Instability detected at a = *([-+0-9.E]+)", log)
    return tag, name, st["status"], why[:2], at[:1], out


def main() -> None:
    jobs = [(harness.TREES["A0"], "A0", n, ov) for n, ov in POINTS.items()]
    jobs += [(harness.TREES["A0"], "A0skip", n, {**ov, "EFTCAMB_skip_stability": "T"}) for n, ov in POINTS.items()]
    jobs += [(R5, "R5", n, ov) for n, ov in POINTS.items()]
    res = {}
    with ThreadPoolExecutor(4) as ex:
        for tag, name, status, why, at, out in ex.map(job, jobs):
            print(f"{tag:8s} {name:8s} {status:12s} {why} {('at a=' + at[0]) if at else ''}")
            res[(tag, name)] = (status, out)
    print("\nIndependent: dtauda argument from cached EFT functions (A0, skip_stability)")
    for name in POINTS:
        status, out = res[("A0skip", name)]
        if status != "ok":
            print(name, "no output", status)
            continue
        d = harness.load(out)
        a = d["eft_a_grid"]
        c, lam = d["eft_EFTc"], d["eft_EFTLambda"]
        om, omp = d["eft_EFTOmegaV"], d["eft_EFTOmegaP"]
        grhom = d["eft_grhom_t"]
        den = 3 * (1 + om + a * omp)
        de_only = (2 * c * a**2 - lam * a**2) / den
        full = (grhom * a**2 + 2 * c * a**2 - lam * a**2) / den
        i = 0
        print(f"{name}: a={a[i]:.1e} c={c[i]:.3e} Lambda={lam[i]:.3e} DE-only arg={de_only[i]:.3e} "
              f"full arg={full[i]:.3e}; n(a) with DE-only arg<0: {(de_only < 0).sum()}/{a.size}, "
              f"min a with DE-only<0: {a[de_only < 0].min() if (de_only < 0).any() else None}")
        j = np.argmin(abs(a - 1e-6))
        print(f"   grhoa2(bridge) at a={a[j]:.2e}: {d['eft_grhoa2'][j]:.4e}  grhom_t*a^2={grhom[j]*a[j]**2:.4e}"
              f"  grhom_t/a^4={grhom[j]/a[j]**4:.4e}")
    st, out = res[("R5", "Om1e-2")]
    if st == "ok":
        d = harness.load(out)
        a = d["eft_a_grid"]
        j = np.argmin(abs(a - 1e-6))
        print(f"R5 bridge grhoa2 at a={a[j]:.2e}: {d['eft_grhoa2'][j]:.4e} vs {d['eft_grhom_t'][j]*a[j]**2:.4e}")
    json.dump({f"{k[0]}_{k[1]}": v[0] for k, v in res.items()}, open(OUT / "status.json", "w"))


if __name__ == "__main__":
    main()

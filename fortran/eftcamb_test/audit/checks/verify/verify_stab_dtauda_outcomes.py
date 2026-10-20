"""Verify t3-stability-dtauda-without-matter (part A): accept/reject matrix.

For GR-limit points of the full-mapping models (OL gamma, OL alpha, K-mouflage, K-mimic) and their nominal test
points, run on A0 (code under audit), A0 with EFTCAMB_skip_stability=T, and R5 (A0 + eft_cache%grhoa2 set in
EFTStabilityComputation). Records status + rejection reason. For runs accepted on R5 but rejected on A0, compare the
spectra with an LCDM reference (EFTflag=0) to see whether the rejected point is physically fine.

Run: python verify_stab_dtauda_outcomes.py
"""

import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REPO = HERE.parents[4]
PARAMS = REPO / "fortran" / "eftcamb_test" / "parameters"
INIS = REPO / "fortran" / "eftcamb_test" / "audit" / "inis"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "stabdtauda" / "outcomes"
R5 = harness.AUDIT_ROOT / "work" / "t3-gr" / "R5"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "1500", "get_transfer": "F",
        "num_redshiftwindows": "0"}

CASES = {
    "OLg_GR": (INIS / "007p5-OL-gamma.ini", {"OLOmegamodel": "0"}),
    "OLg_1e-6": (INIS / "007p5-OL-gamma.ini", {"OLOmega0": "1e-6"}),
    "OLg_nom": (INIS / "007p5-OL-gamma.ini", {}),
    "OLa_GR": (INIS / "007p6-OL-alpha.ini", {"OLmassPmodel": "0", "OLkineticitymodel": "0"}),
    "OLa_1e-6": (INIS / "007p6-OL-alpha.ini", {"OLmass0": "1e-6"}),
    "OLa_nom": (INIS / "007p6-OL-alpha.ini", {}),
    "Kmou_GR": (PARAMS / "5_Kmouflage_1.ini", {"eps2_0": "0"}),
    "Kmou_-1e-6": (PARAMS / "5_Kmouflage_1.ini", {"eps2_0": "-1e-6"}),
    "Kmou_nom": (PARAMS / "5_Kmouflage_1.ini", {}),
    "Kmim_GR": (PARAMS / "5_Kmimic_1.ini", {"eps2_0": "0"}),
    "Kmim_1e-6": (PARAMS / "5_Kmimic_1.ini", {"eps2_0": "1e-6"}),
    "Kmim_nom": (PARAMS / "5_Kmimic_1.ini", {}),
}
TREES = {"A0": "A0", "A0skip": "A0", "R5": str(R5)}
REASON = re.compile(r"(Model dtauda is Nan|Physical instability[^\n]*|Mathematical instability[^\n]*|"
                    r"Model has Nan[^\n]*|Error in Fortran[^\n]*)")
AT = re.compile(r"Instability detected at a = *([-+0-9.E]+)")


def job(item):
    name, (tree, ini, ov) = item
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, ini, out, overrides=ov, threads=1, timeout=540)
    log = out.with_suffix(".log").read_text(errors="replace")
    return name, st["status"], REASON.findall(log)[:1], AT.findall(log)[:1]


def main() -> None:
    jobs = {}
    for c, (ini, ov) in CASES.items():
        for t, tree in TREES.items():
            o = {**FAST, **ov}
            if t == "A0skip":
                o["EFTCAMB_skip_stability"] = "T"
            jobs[f"{t}__{c}"] = (tree, ini, o)
    # LCDM references (EFTflag=0) with and without massive neutrinos (K-mouflage/K-mimic inis set omnuh2=0)
    jobs["REF__lcdm"] = ("A0", PARAMS / "base_params.ini", {**FAST, "EFTflag": "0"})
    jobs["REF__lcdm_nu0"] = ("A0", PARAMS / "base_params.ini", {**FAST, "EFTflag": "0", "omnuh2": "0.0"})
    with ThreadPoolExecutor(4) as ex:
        res = {n: (s, w, a) for n, s, w, a in ex.map(job, jobs.items())}
    table = {}
    for n, (s, w, a) in res.items():
        line = f"{n:22s} {s:12s} {w} {('at a=' + a[0]) if a else ''}"
        if s == "ok" and not n.startswith("REF"):
            case = n.split("__")[1]
            ref = "REF__lcdm_nu0" if case.startswith("K") else "REF__lcdm"
            if res[ref][0] == "ok":
                cmp = harness.compare(OUT / f"{n}.npz", OUT / f"{ref.split('__')[0]}__{ref.split('__')[1]}.npz",
                                      lmax=1500)
                line += f"  vs LCDM: worst_cl={cmp['worst_cl']:.3e}"
        print(line)
        table[n] = {"status": s, "reason": w, "at": a}
    # A0 vs R5 identity on accepted nominal points (the fix must not change spectra of accepted models)
    for c in CASES:
        a0, r5 = OUT / f"A0__{c}.npz", OUT / f"R5__{c}.npz"
        if res[f"A0__{c}"][0] == "ok" and res[f"R5__{c}"][0] == "ok":
            cmp = harness.compare(a0, r5, lmax=1500)
            print(f"A0 vs R5 {c:12s} bit_identical={cmp['bit_identical']} worst_cl={cmp['worst_cl']:.3e}")
    (OUT / "summary.json").write_text(json.dumps(table, indent=1))


if __name__ == "__main__":
    main()

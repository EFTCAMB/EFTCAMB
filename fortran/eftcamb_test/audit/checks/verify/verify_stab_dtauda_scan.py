"""Verify t3-stability-dtauda-without-matter (part C): does the missing grhoa2 ever flip accept/reject?

Scan parameter points of models whose compute_dtauda depends on eft_cache%grhoa2 (OL gamma / OL alpha full-mapping
models, Horava designer) on A0 and on R5 (A0 + eft_cache%grhoa2 = grhom_t a^2 in EFTStabilityComputation). A "flip"
is a point rejected on A0 with 'Model dtauda is Nan' and accepted on R5. Flipped points are also run on A0 with
EFTCAMB_skip_stability=T to check that the R5 spectra equal the unchecked A0 spectra (the fix only touches the check).

Run: python verify_stab_dtauda_scan.py [group ...]   (groups: horava, olg, ola)
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
OUT = harness.AUDIT_ROOT / "work" / "verify" / "stabdtauda" / "scan"
R5 = harness.AUDIT_ROOT / "work" / "t3-gr" / "R5"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "600", "get_transfer": "F",
        "num_redshiftwindows": "0"}
REASON = re.compile(r"(Model dtauda is Nan|Physical instability[^\n]*|Mathematical instability[^\n]*|"
                    r"Model has Nan[^\n]*|Recfast: [^\n]*|theory unstable|theory stable)")


def groups() -> dict:
    g = {"horava": {}, "olg": {}, "ola": {}}
    hor = PARAMS / "5_Horava_eta_1.ini"
    for xi, lam, eta in [(2.0, 0.0, 0.0), (2.5, 0.0, 0.0), (3.0, 0.0, 0.0), (5.0, 0.0, 0.0), (0.0, -0.5, 0.0),
                         (0.0, -0.6, 0.0), (0.0, 0.0, -1.0), (0.0, 0.0, 1.0)]:
        ov = {"Horava_xi": repr(xi), "Horava_lambda": repr(lam), "Horava_eta": repr(eta)}
        g["horava"][f"hor_xi{xi}_l{lam}_e{eta}"] = (hor, ov)
        # same with the physical ghost and gradient checks on
        g["horava"][f"horphys_xi{xi}_l{lam}_e{eta}"] = (hor, {**ov, "EFT_ghost_stability": "T",
                                                              "EFT_gradient_stability": "T"})
    olg = INIS / "007p5-OL-gamma.ini"
    for om in ("-0.05", "-0.01", "0.05"):
        for g1 in ("0", "0.1", "1.0"):
            for g2 in ("0", "0.1", "-0.1"):
                ov = {"OLOmega0": om}
                if g1 != "0":
                    ov.update({"OLGamma1model": "1", "OLGamma10": g1})
                if g2 != "0":
                    ov.update({"OLGamma2model": "1", "OLGamma20": g2})
                g["olg"][f"olg_om{om}_g1{g1}_g2{g2}"] = (olg, ov)
    ola = INIS / "007p6-OL-alpha.ini"
    for m0 in ("-0.05", "-0.01", "0.05", "0.2"):
        for k0 in ("1.0", "10.0"):
            for b0 in ("0", "0.1", "-0.1"):
                ov = {"OLmass0": m0, "OLkineticity0": k0}
                if b0 != "0":
                    ov.update({"OLbraidingmodel": "1", "OLbraiding0": b0})
                g["ola"][f"ola_m{m0}_k{k0}_b{b0}"] = (ola, ov)
    return g


def job(item):
    name, (tree, ini, ov) = item
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, ini, out, overrides={**FAST, **ov}, threads=1, timeout=300)
    log = out.with_suffix(".log").read_text(errors="replace")
    rs = [r.strip() for r in dict.fromkeys(REASON.findall(log)) if r not in ("theory unstable", "theory stable")]
    return name, st["status"], rs[:2]


def main() -> None:
    sel = sys.argv[1:] or ["horava", "olg", "ola"]
    allg = groups()
    cases = {k: v for s in sel for k, v in allg[s].items()}
    jobs = {f"{t}__{c}": (tree, ini, ov) for c, (ini, ov) in cases.items()
            for t, tree in (("A0", "A0"), ("R5", str(R5)))}
    with ThreadPoolExecutor(4) as ex:
        res = {n: (s, r) for n, s, r in ex.map(job, jobs.items())}
    flips = []
    for c in cases:
        a, r = res[f"A0__{c}"], res[f"R5__{c}"]
        flip = a[0] != "ok" and "Model dtauda is Nan" in a[1] and r[0] == "ok"
        print(f"{c:34s} A0: {a[0]:11s} {str(a[1])[:60]:60s} R5: {r[0]:11s} {str(r[1])[:70]}{'  <== FLIP' if flip else ''}")
        if flip:
            flips.append(c)
    skip_jobs = {f"A0skip__{c}": ("A0", cases[c][0], {**cases[c][1], "EFTCAMB_skip_stability": "T"}) for c in flips}
    with ThreadPoolExecutor(4) as ex:
        list(ex.map(job, skip_jobs.items()))
    for c in flips:
        cmp = harness.compare(OUT / f"R5__{c}.npz", OUT / f"A0skip__{c}.npz", lmax=600)
        print(f"flip {c}: R5 vs A0skip worst_cl={cmp['worst_cl']:.3e}")
    (OUT / f"scan_{'_'.join(sel)}.json").write_text(json.dumps({"res": res, "flips": flips}, indent=1))
    print(f"flips: {len(flips)} / {len(cases)}")


if __name__ == "__main__":
    main()

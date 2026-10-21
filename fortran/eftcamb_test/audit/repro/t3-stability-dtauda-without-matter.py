"""Repro (t3-gr): the stability check evaluates dtauda with zero matter density, rejecting GR-limit models.

09_EFTCAMB_stability.f90 EFTTestStability calls model%compute_dtauda(a, params_cache, eft_cache) right after
EFTStabilityComputation, which fills grhom_t, grhob_t, ... but never eft_cache%grhoa2 (it stays 0 from
eft_cache%initialize()). Every compute_dtauda implementation adds the matter density through eft_cache%grhoa2
(e.g. EFTCAMBFullModelComputeDtauda: 1/sqrt((grhoa2 + 2 c a^2 - Lambda a^2)/3/(1+Omega+a Omega'))), so the check
tests the sign of the dark-energy part alone. Where 2c - Lambda < 0 (here: the GR point of the Omega-Lambda gamma
model, whose background solver leaves c ~ -2.6e-4 at a = 1e-9, while Lambda a^2 ~ 1e-25) the result is NaN and the
model is rejected ("Model dtauda is Nan" at a = 1e-10, "EFTCAMB: theory unstable") although H^2 > 0.

Model: OL gamma (AltParEFTmodel=2, ini fortran/eftcamb_test/audit/inis/007p5-OL-gamma.ini) approaching GR,
OLOmega0 = 1e-2, 1e-4, 1e-6, 0 (OLOmegamodel=0). Control: same with EFTCAMB_skip_stability=T.
Optional: a tree with eft_cache%grhoa2 = eft_cache%grhom_t*a**2 added to EFTStabilityComputation (pass its path).
Expected: all points accepted (they converge to LCDM). Observed: printed.
Run: /Users/marco/Software/miniforge3/bin/python t3-stability-dtauda-without-matter.py [FIXED_TREE]
"""

import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

INI = REPO / "fortran" / "eftcamb_test" / "audit" / "inis" / "007p5-OL-gamma.ini"
OUT = harness.AUDIT_ROOT / "work" / "t3-gr" / "repro_stab_dtauda"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "100", "get_transfer": "F",
        "num_redshiftwindows": "0"}
POINTS = {"Om1e-2": {"OLOmega0": "1e-2"}, "Om1e-4": {"OLOmega0": "1e-4"}, "Om1e-6": {"OLOmega0": "1e-6"},
          "Om0": {"OLOmegamodel": "0"}}


def job(args):
    tree, tag, name, ov = args
    out = OUT / f"{tag}_{name}.npz"
    st = harness.run_model(tree, INI, out, overrides={**FAST, **ov}, threads=1, timeout=900)
    log = out.with_suffix(".log").read_text(errors="replace")
    why = re.findall(r"(Model dtauda is Nan|Physical instability[^\n]*|Mathematical instability[^\n]*)", log)
    at = re.findall(r"Instability detected at a = *([-+0-9.E]+)", log)
    return tag, name, st["status"], why[:1], at[:1]


def main() -> None:
    jobs = [("A0", "A0", n, ov) for n, ov in POINTS.items()]
    jobs += [("A0", "A0_skipstab", n, {**ov, "EFTCAMB_skip_stability": "T"}) for n, ov in POINTS.items()]
    if len(sys.argv) > 1:
        jobs += [(sys.argv[1], "fixed", n, ov) for n, ov in POINTS.items()]
    with ThreadPoolExecutor(4) as ex:
        for tag, name, status, why, at in ex.map(job, jobs):
            print(f"{tag:12s} {name:8s} {status:12s} {why} {('at a=' + at[0]) if at else ''}")
    print("expected: status ok for every point (the model tends to LCDM); A0 rejects OLOmega0<=1e-6 and the GR point")


if __name__ == "__main__":
    main()

"""Repro (t3-gr): the 3+1 wDE designer (EFTflag=3, DesignerEFTmodel=3) always fails with "background solver failed".

07f_designer_models/007p8_3p1_wDE_designer.f90: EFTCAMBDesigner3p1InitBackground sets success = .true. and passes it
to solve_designer_equations, whose dummy is `logical, intent(out) :: success` and is only assigned on the DLSODA
error path (success = .false.). On the success path the intent(out) dummy is undefined on return; the A0 release
build returns .false. for every parameter point, including the GR point (sigma1 = rho1 = 1, w = -1), so no 3+1
model can be run. (Same code in the pre-rebase B1, where it happened to return .true.)
Fix: set `success = .true.` at the top of EFTCAMBDesigner3p1SolveDesignerEquations.
With the fix (tree work/t3-gr/R7, optional argument) the GR point with EFTCAMB_GR_threshold=-1 and
EFTCAMB_skip_stability=T reproduces EFTflag=0 (IntTolBoost=2) to 1.2e-8 in C_l (with the stability check on, the
GR point is then still rejected as "beyond GLPV" because EFTGamma3V = sigma1/3 - rho1 + 2/3 = -1.1e-16 by rounding).
Expected: status ok. Observed: printed.
Run: /Users/marco/Software/miniforge3/bin/python t3-3p1-success-unset.py [FIXED_TREE]
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

INI = REPO / "fortran" / "eftcamb_test" / "audit" / "inis" / "007p8-3p1.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
OUT = harness.AUDIT_ROOT / "work" / "t3-gr" / "repro_3p1"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "200", "get_transfer": "F",
        "num_redshiftwindows": "0"}
GRPT = {"EFTwDE": "0", "Designer3p1modelsigma1": "1", "Designer3p1sigma10": "1.0", "Designer3p1modelsigma1_ODE": "0",
        "Designer3p1modelrho1N": "0", "Designer3p1modelrho2": "0"}
FORCED = {"EFTCAMB_GR_threshold": "-1", "EFTCAMB_skip_stability": "T"}


def job(args):
    tree, name, ini, ov, ps = args
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, ini, out, overrides={**FAST, **ov}, python_set=ps, threads=1, timeout=600)
    msg = [ln.strip() for ln in out.with_suffix(".log").read_text(errors="replace").splitlines()
           if "background solver failed" in ln or "theory unstable" in ln or "beyond GLPV" in ln][-1:]
    return name, st, out, msg


def main() -> None:
    jobs = [("A0", "A0_3p1_ini", INI, {}, {}),
            ("A0", "A0_3p1_GRpoint", INI, GRPT, {}),
            ("A0", "A0_3p1_GRpoint_forced", INI, {**GRPT, **FORCED}, {}),
            ("A0", "A0_LCDM", GR, {}, {"Accuracy.IntTolBoost": 2.0})]
    if len(sys.argv) > 1:
        jobs += [(sys.argv[1], "fixed_3p1_GRpoint_forced", INI, {**GRPT, **FORCED}, {}),
                 (sys.argv[1], "fixed_3p1_GRpoint", INI, GRPT, {})]
    with ThreadPoolExecutor(3) as ex:
        res = {n: (st, out, msg) for n, st, out, msg in ex.map(job, jobs)}
    ref = res["A0_LCDM"][1]
    for n, (st, out, msg) in res.items():
        if n == "A0_LCDM":
            continue
        line = f"{n:26s} {st['status']:12s} {msg}"
        if st["status"] == "ok":
            line += f"  max|dCl/Cl| vs EFTflag=0: {harness.compare(out, ref)['worst_cl']:.1e}"
        print(line)
    print("expected: status ok (and C_l = EFTflag=0 at the GR point); observed: 'background solver failed' on A0")


if __name__ == "__main__":
    main()

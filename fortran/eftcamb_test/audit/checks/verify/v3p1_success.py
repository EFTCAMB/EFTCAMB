"""Adversarial verifier for t3-3p1-success-unset (F-006).

Stage 'runs': A0 vs fixed tree (FIX = A0 + `success = .true.` at top of EFTCAMBDesigner3p1SolveDesignerEquations),
for the 3+1 test ini, its GR point and the GR point with forced EFT path; plus EFTflag=0 reference.
Stage 'build': build the fixed tree.
Run with AUDIT_ROOT set to the live scratchpad audit root.
"""

import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

INI = HERE / "v3p1_success.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "v3p1"
FIX = harness.AUDIT_ROOT / "work" / "verify" / "v3p1" / "FIX"
SRC = "fortran/eftcamb/07f_designer_models/007p8_3p1_wDE_designer.f90"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "200", "get_transfer": "F",
        "num_redshiftwindows": "0"}
GRPT = {"EFTwDE": "0", "Designer3p1modelsigma1": "1", "Designer3p1sigma10": "1.0", "Designer3p1modelsigma1_ODE": "0",
        "Designer3p1modelrho1N": "0", "Designer3p1modelrho2": "0"}
FORCED = {"EFTCAMB_GR_threshold": "-1", "EFTCAMB_skip_stability": "T"}


def build():
    if not FIX.exists():
        import shutil
        shutil.copytree(harness.TREES["A0"], FIX, symlinks=True)
    f = FIX / SRC
    txt = f.read_text()
    old = "    ! 1) Cosmological parameters:\n    grhom0 = params_cache%grhob + params_cache%grhoc\n    grhor0 = 3._dl*params_cache%h0_Mpc**2*(params_cache%omegag + params_cache%omegar)\n    grhov0 = params_cache%grhov\n\n    ! 2) Initial values\n    y(1) = 1._dl"
    if "success = .true. ! VERIFYFIX" not in txt:
        assert txt.count(old) == 1
        txt = txt.replace(old, "    success = .true. ! VERIFYFIX\n" + old)
        f.write_text(txt)
    print(harness.build_tree(FIX, FIX))


def job(args):
    tree, name, ini, ov, ps = args
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, ini, out, overrides={**FAST, **ov}, python_set=ps, threads=1, timeout=600)
    log = out.with_suffix(".log")
    msg = []
    if log.exists():
        msg = [ln.strip() for ln in log.read_text(errors="replace").splitlines()
               if any(s in ln for s in ("solver failed", "unstable", "GLPV", "DLSODA", "Error", "error"))][-2:]
    return name, st, out, msg


def runs():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [("A0", "A0_ini", INI, {}, {}),
            ("A0", "A0_GR", INI, GRPT, {}),
            ("A0", "A0_GRforced", INI, {**GRPT, **FORCED}, {}),
            ("B1", "B1_ini", INI, {}, {}),
            ("A0", "A0_LCDM", GR, {}, {"Accuracy.IntTolBoost": 2.0})]
    if (FIX / "camb" / "camblib.so").exists():
        jobs += [(str(FIX), "FIX_ini", INI, {}, {}),
                 (str(FIX), "FIX_GR", INI, GRPT, {}),
                 (str(FIX), "FIX_GRforced", INI, {**GRPT, **FORCED}, {})]
    with ThreadPoolExecutor(3) as ex:
        res = {n: (st, out, msg) for n, st, out, msg in ex.map(job, jobs)}
    ref = res["A0_LCDM"][1]
    for n, (st, out, msg) in res.items():
        line = f"{n:14s} {st['status']:12s} {msg}"
        if st["status"] == "ok" and n != "A0_LCDM":
            line += f"  worst dCl vs EFTflag=0: {harness.compare(out, ref)['worst_cl']:.2e}"
        print(line)


if __name__ == "__main__":
    assert os.environ.get("AUDIT_ROOT"), "set AUDIT_ROOT"
    {"build": build, "runs": runs}[sys.argv[1]]()

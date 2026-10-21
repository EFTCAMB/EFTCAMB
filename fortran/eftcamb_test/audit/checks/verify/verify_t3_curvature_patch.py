"""Attribution test for t3-eft-curvature-ignored: a DIAGNOSTIC patch (not a proposed fix) of a copy of A0.

In derivs (equations.f90, EFT branch, after all EFT background functions are computed with the flat H), restore the
curvature in the conformal Hubble rate, in adotdota (CAMB's GR expression) and in cothxor. If the stated root cause is
right, PEFT-zero with EFT perturbations never switched on (turn_on_time=1.1) must then agree with EFTflag=0 at
omk/=0 at the 1e-5 level. With the EFT perturbations active (turn_on=0.01) the EFT Einstein equations still lack the
curvature terms, so a residual is expected there.

Usage: python verify_t3_curvature_patch.py build | run
"""

import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import harness  # noqa: E402
import verify_t3_curvature as V  # noqa: E402

DEST = harness.AUDIT_ROOT / "work" / "verify" / "curvdiag_tree"

OLD = """        ! call CP%EFTCAMB%model%EFTCAMBModelComputeQSMuSigma( a, k, CP%eft_par_cache , EV%eft_cache ) minsup
        !
        cothxor = 1._dl/tau
"""
NEW = """        ! call CP%EFTCAMB%model%EFTCAMBModelComputeQSMuSigma( a, k, CP%eft_par_cache , EV%eft_cache ) minsup
        !
        cothxor = 1._dl/tau
        ! AUDIT DIAGNOSTIC: curvature in H, adotdota and cothxor
        if (.not. State%flat) then
            adotoa = sqrt(adotoa**2 + State%grhok/3._dl)
            adotdota = (adotoa*adotoa - gpres_noDE - EV%eft_cache%gpiv_t)/2._dl
            cothxor = 1._dl/State%tanfunc(tau/State%curvature_radius)/State%curvature_radius
        end if
"""


OLD2 = """            if (CP%EFTCAMB%EFTflag == 0 .or. .not. EV%EFTCAMBactive) then
                if (CP%EFTCAMB%EFTflag /= 0) then
                    ISW = -2._dl*sigma*(EV%eft_cache%Hdot - 2._dl*adotoa**2) - 2._dl*adotoa*etak + dgq - diff_rhopi/k
"""
NEW2 = OLD2.replace("if (CP%EFTCAMB%EFTflag /= 0) then", "if (CP%EFTCAMB%EFTflag /= 0 .and. State%flat) then ! AUDIT DIAG2")


def build(stage: str = "1") -> None:
    """stage 1: H/adotdota/cothxor patch only; stage 2: also CAMB's curved ISW (2 phidot) when EFT is inactive."""
    if not DEST.exists():
        shutil.copytree(harness.tree_path("A0"), DEST, symlinks=True)
    f = DEST / "fortran" / "equations.f90"
    s = f.read_text()
    if "AUDIT DIAGNOSTIC" not in s:
        assert s.count(OLD) == 1
        s = s.replace(OLD, NEW)
    if stage == "2" and "AUDIT DIAG2" not in s:
        assert s.count(OLD2) == 1
        s = s.replace(OLD2, NEW2)
    f.write_text(s)
    for _ in range(4):  # forutils' Makefile races under make -j4 on a fresh build dir: retry
        r = harness.build_tree("A0", DEST)
        if "camblib.so ../camb/" in Path(r["log"]).read_text():
            break
    print(r)


def run() -> None:
    V.OUT.mkdir(parents=True, exist_ok=True)
    tree = str(DEST)
    jobs = [
        ("D_PEFT_on1.1_omk0.05", tree, 0.05, {**V.PEFT_ZERO, "EFTCAMB_turn_on_time": 1.1}),
        ("D_PEFT_on1.1_omk-0.05", tree, -0.05, {**V.PEFT_ZERO, "EFTCAMB_turn_on_time": 1.1}),
        ("D_PEFT_omk0.05", tree, 0.05, V.PEFT_ZERO),
        ("D_PEFT_omk0.0", tree, 0.0, V.PEFT_ZERO),
        ("D_GR_omk0.05", tree, 0.05, None),
        ("A0_GR_omk-0.05", "A0", -0.05, None),
    ]
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(4) as ex:
        st = list(ex.map(lambda j: V.run(*j), jobs))
    for s in st:
        if s["rc"] != 0:
            print("FAILED", s)
    rows = [
        ("patched PEFT turn_on=1.1 vs A0 GR, omk=+0.05", "D_PEFT_on1.1_omk0.05", "A0_GR_omk0.05"),
        ("patched PEFT turn_on=1.1 vs A0 GR, omk=-0.05", "D_PEFT_on1.1_omk-0.05", "A0_GR_omk-0.05"),
        ("patched PEFT turn_on=0.01 vs A0 GR, omk=+0.05", "D_PEFT_omk0.05", "A0_GR_omk0.05"),
        ("patched PEFT vs A0 GR, omk=0 (patch inert)", "D_PEFT_omk0.0", "A0_GR_omk0.0"),
        ("patched GR vs A0 GR, omk=0.05 (patch inert)", "D_GR_omk0.05", "A0_GR_omk0.05"),
        ("unpatched PEFT turn_on=1.1 vs A0 GR, omk=+0.05", "A0_PEFT_on1.1_omk0.05", "A0_GR_omk0.05"),
    ]
    keys = ("TT", "TT_l>=30", "EE", "TE", "pk", "sigma8", "Hbg")
    print(f"{'case':50s} " + " ".join(f"{k:>9s}" for k in keys))
    for lab, a, b in rows:
        r = V.summary(a, b)
        print(f"{lab:50s} " + " ".join(f"{r[k]:9.2e}" for k in keys))


if __name__ == "__main__":
    if sys.argv[1] == "build":
        build(sys.argv[2] if len(sys.argv) > 2 else "1")
    else:
        run()

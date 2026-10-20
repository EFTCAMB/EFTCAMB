"""Build a modified copy of a frozen audit tree (t3-gr): copy tree, apply exact string patches, build in own dir.

Usage: python build_variant.py NAME [--base A0] [--fflags "..."] PATCH [PATCH ...]
PATCH names are keys of PATCHES below. The copy goes to $WORK/NAME, the build dir to $WORK/NAME_build.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import AUDIT, WORK  # noqa: E402

PY = "/Users/marco/Software/miniforge3/bin/python"
COMMON = "-MMD -cpp -ffree-line-length-none -fallow-argument-mismatch -fmax-errors=4 -fopenmp -lblas -llapack"

POLTER_NEW = """            pigddot = 2._dl/5._dl*k*qgdot - 3._dl/5._dl*k*octgdot - dopacity*(pig - polter) - &
                opacity*(pigdot - polterdot) + 8._dl/15._dl*k*sigmadot
            E2ddot = -dopacity*(E(2) - polter) - opacity*(Edot(2) - polterdot) - k/3._dl*Edot(3)
            polterddot = pigddot/10._dl + 3._dl/5._dl*E2ddot
"""
POLTER_OLD = """            polterddot = -2.0d0/25.0d0*adotoa*dgq/(k*EV%Kf(1)) - 4.0d0/75.0d0*adotoa* &
                k*sigma - 4.0d0/75.0d0*dgpi - 2.0d0/75.0d0*dgrho/EV%Kf(1) - 3.0d0/ &
                50.0d0*k*octgdot*EV%Kf(2) + (1.0d0/25.0d0)*k*qgdot - 1.0d0/5.0d0 &
                *k*EV%Kf(2)*Edot(3) + (-1.0d0/10.0d0*pig + (7.0d0/10.0d0)* &
                polter - 3.0d0/5.0d0*E(2))*dopacity + (-1.0d0/10.0d0*pigdot &
                + (7.0d0/10.0d0)*polterdot - 3.0d0/5.0d0*Edot(2))*opacity
"""
# curvature-consistent version of the EFTCAMB photon-equation polterddot (Kf(2) on the octupole term, as in pigdot)
POLTER_KF2 = """            pigddot = 2._dl/5._dl*k*qgdot - 3._dl/5._dl*k*EV%Kf(2)*octgdot - dopacity*(pig - polter) - &
                opacity*(pigdot - polterdot) + 8._dl/15._dl*k*sigmadot
            E2ddot = -dopacity*(E(2) - polter) - opacity*(Edot(2) - polterdot) - k/3._dl*Edot(3)
            polterddot = pigddot/10._dl + 3._dl/5._dl*E2ddot
"""

# full fix: Kf(2) on both the octupole term (pigdot) and the E3 term, as in upstream CAMB's polterddot and in
# camb/symbolic.py's E hierarchy (algebraically identical to the original CAMB formula in GR)
POLTER_KF22 = """            pigddot = 2._dl/5._dl*k*qgdot - 3._dl/5._dl*k*EV%Kf(2)*octgdot - dopacity*(pig - polter) - &
                opacity*(pigdot - polterdot) + 8._dl/15._dl*k*sigmadot
            E2ddot = -dopacity*(E(2) - polter) - opacity*(Edot(2) - polterdot) - k/3._dl*EV%Kf(2)*Edot(3)
            polterddot = pigddot/10._dl + 3._dl/5._dl*E2ddot
"""

# name -> list of (file, old, new)
PATCHES: dict[str, list[tuple[str, str, str]]] = {
    "maxz50": [(
        "fortran/reionization.f90",
        "        real(dl) :: max_redshift = 1000._dl\n",
        "        real(dl) :: max_redshift = 50._dl\n",
    )],
    "polter_orig": [("fortran/equations.f90", POLTER_NEW, POLTER_OLD)],
    "polter_kf2": [("fortran/equations.f90", POLTER_NEW, POLTER_KF2)],
    "polter_kf22": [("fortran/equations.f90", POLTER_NEW, POLTER_KF22)],
    "weyl_orig": [(
        "fortran/equations.f90",
        "                EV%OutputTransfer(Transfer_Weyl) = -(dgrho + 3*dgq*adotoa/k)/(EV%Kf(1)*2) - dgpi/2\n",
        "                EV%OutputTransfer(Transfer_Weyl) = k2*phi\n",
    )],
    # stability check: fill grhoa2 (= 8 pi G rho_matter a^4) so compute_dtauda sees the matter density
    "stab_grhoa2": [(
        "fortran/eftcamb/09_EFTCAMB_stability.f90",
        "    eft_cache%gpresm_t = (+eft_cache%grhor_t + eft_cache%grhog_t)/3._dl + eft_cache%gpinu_tot\n"
        "    ! compute the other things:\n",
        "    eft_cache%gpresm_t = (+eft_cache%grhor_t + eft_cache%grhog_t)/3._dl + eft_cache%gpinu_tot\n"
        "    eft_cache%grhoa2 = eft_cache%grhom_t*a**2\n"
        "    ! compute the other things:\n",
    )],
    # python bridge (get_eft_functions): grhoa2 = grhom_t*a^2, not grhom_t/a^4
    "py_grhoa2": [(
        "fortran/eftcamb/11_EFTCAMB_python.f90",
        "                cache%grhoa2 = cache%grhom_t/a**4\n",
        "                cache%grhoa2 = cache%grhom_t*a**2\n",
    )],
    # number-count RSD source in the GR branch with EFTflag/=0: grho/gpres lack the (EFT) dark energy; use the
    # Friedmann equations (gpres_tot + grho_tot/3 = -2 Hdot) instead
    "counts_rsd": [(
        "fortran/equations.f90",
        "                        counts_redshift_source = ((4.d0*adotoa**2 + gpres + grho/3.d0)/k*w%wing2(j) + &\n",
        "                        counts_redshift_source = ((4.d0*adotoa**2 + merge(gpres + grho/3.d0, &\n"
        "                            -2._dl*EV%eft_cache%Hdot, CP%EFTCAMB%EFTflag == 0))/k*w%wing2(j) + &\n",
    )],
    # 3+1 designer: SolveDesignerEquations has intent(out) success that is never set on the success path
    "p31_success": [(
        "fortran/eftcamb/07f_designer_models/007p8_3p1_wDE_designer.f90",
        "    grhom0 = params_cache%grhob + params_cache%grhoc\n",
        "    success = .true.\n    grhom0 = params_cache%grhob + params_cache%grhoc\n",
    )],
    # pi-dot initial condition: d/da of (D1 + k^2 D2), not (D1 + k^2 D1)
    "pid_d2": [(
        "fortran/eftcamb/09_EFTCAMB_IC.f90",
        "    EFTpiDfunction = eft_cache%EFTpiD1 + k*k*eft_cache%EFTpiD1\n",
        "    EFTpiDfunction = eft_cache%EFTpiD1 + k*k*eft_cache%EFTpiD2\n",
    )],
    # get_evolution after a transfer-only run: GetNumEqns sets up the scalar hierarchy only if CP%WantScalars
    "getevol_scalars": [
        ("fortran/eftcamb/11_EFTCAMB_python.f90",
         "    real(dl) :: taustart\n    integer  :: q_ix\n    type(EvolutionVars) :: Ev\n",
         "    real(dl) :: taustart\n    integer  :: q_ix\n    type(EvolutionVars) :: Ev\n    logical  :: want_scalars\n"),
        ("fortran/eftcamb/11_EFTCAMB_python.f90",
         "    ! do the parallel loop over the requested k modes.\n",
         "    want_scalars = this%CP%WantScalars\n    this%CP%WantScalars = .true.\n"
         "    ! do the parallel loop over the requested k modes.\n"),
        ("fortran/eftcamb/11_EFTCAMB_python.f90",
         "    !$OMP END PARALLEL DO\n\n    end subroutine EFTCAMB_GetEvolution\n",
         "    !$OMP END PARALLEL DO\n    this%CP%WantScalars = want_scalars\n\n    end subroutine EFTCAMB_GetEvolution\n"),
    ],
    # upstream F-002 workaround needed for any -fcheck build
    "f002": [(
        "fortran/results.f90",
        "    if (allocated(this%tb) .and. this%nthermo /= size(this%tb)) then",
        "    if (allocated(this%tb)) then\n    if (this%nthermo /= size(this%tb)) then",
    )],
}


def apply(tree: Path, name: str) -> None:
    for rel, old, new in PATCHES[name]:
        f = tree / rel
        txt = f.read_text()
        n = txt.count(old)
        if n != 1:
            raise SystemExit(f"patch {name}: expected 1 match in {rel}, found {n}")
        f.write_text(txt.replace(old, new))
        print(f"applied {name} to {rel}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("name")
    p.add_argument("patches", nargs="*")
    p.add_argument("--base", default="A0")
    p.add_argument("--fflags", default=None, help="full FFLAGS (COMMON is appended)")
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    src = AUDIT / "trees" / a.base
    dst = WORK / a.name
    bld = WORK / f"{a.name}_build"
    if dst.exists():
        if not a.force:
            raise SystemExit(f"{dst} exists (use --force)")
        shutil.rmtree(dst)
        shutil.rmtree(bld, ignore_errors=True)
    shutil.copytree(src, dst, symlinks=True)
    (dst / "camb" / "camblib.so").unlink(missing_ok=True)
    for name in a.patches:
        apply(dst, name)
    info = json.loads((src / "AUDIT_TREE.json").read_text())
    info.update({"variant_of": a.base, "patches": a.patches, "fflags": a.fflags})
    (dst / "AUDIT_TREE.json").write_text(json.dumps(info, indent=2))
    # one build at a time machine-wide (harness lock), make -j4, nice 10; build dir = <dst>_build
    from common import harness

    fflags = None if a.fflags is None else f"{a.fflags} {COMMON}"
    res = harness.build_tree(dst, dst, fflags=fflags)
    ok, rc, log = res["ok"], res["returncode"], res["log"]
    print(f"{a.name}: {'OK' if ok else f'FAILED rc={rc}'} (log {log})")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

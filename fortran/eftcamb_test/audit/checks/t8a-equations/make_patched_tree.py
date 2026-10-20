"""Copy the frozen A0 tree, apply runtime-switchable audit patches to fortran/equations.f90 and build it.

All patches are guarded by logical flags read from environment variables in GaugeInterface_Init, so with no
AUDIT_* variable set the patched build computes exactly what A0 computes (checked by the caller).

Flags (set to 1 to enable):
  AUDIT_FIX_PIDOTDOT  output(): call derivs twice so that eft_cache%pidotdot (read from ayprime before it is
                      computed) is the current pi'' when the Einstein factors (EFTeomM, EFTeomNdot) are evaluated
                      for the line-of-sight sources.
  AUDIT_FIX_DPMATTER  derivs, RSA branch of the delta_phi equations: dp_matter = dgpnu/a2 (as in the non-RSA branch).
  AUDIT_FIX_TENSOR    derivst: set grhog_t, grhor_t before the EFT block that uses them.
  AUDIT_POLTER_GR     derivs output: for EFTflag=0 use the original CAMB polterddot expression (with Kf factors).
  AUDIT_FIX_METRIC_H  derivs: set ayprime(w_ix+2)=0 while EFT is inactive (evolve_delta_phi + evolve_metric_h).
  AUDIT_DGP_EXACT     MassiveNuVars: exact massive-neutrino pressure perturbation instead of p_nu*delta_nu.
  AUDIT_FIX_WDE       derivs: define w_dark_energy_t = -1 for EFTflag /= 0 (it is used as w*grhov_t with grhov_t=0).

Usage: python make_patched_tree.py NAME [release|snan]
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
    "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
)
WORK = ROOT / "work" / "t8a-equations"
PY = "/Users/marco/Software/miniforge3/bin/python"
COMMON = "-MMD -cpp -ffree-line-length-none -fallow-argument-mismatch -fmax-errors=4 -fopenmp -lblas -llapack"
SNAN = "-g -O0 -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-999999 -finit-logical=true"

FLAGS = ["pidotdot", "dpmatter", "tensor", "polter_gr", "metric_h", "dgp_exact", "wde"]

PATCHES = [
    # module flags
    (
        "    real(dl) nu_tau_nonrelativistic_physical(max_nu), nu_tau_nonrelativistic(max_nu), nu_tau_massive(max_nu)\n",
        "    real(dl) nu_tau_nonrelativistic_physical(max_nu), nu_tau_nonrelativistic(max_nu), nu_tau_massive(max_nu)\n"
        "    logical :: audit_fix_pidotdot = .false., audit_fix_dpmatter = .false., audit_fix_tensor = .false., &\n"
        "        audit_polter_gr = .false., audit_fix_metric_h = .false., audit_dgp_exact = .false., &\n"
        "        audit_fix_wde = .false.\n",
    ),
    (
        "    epsw = 100/State%tau0\n",
        "    epsw = 100/State%tau0\n"
        "    audit_fix_pidotdot = audit_env('AUDIT_FIX_PIDOTDOT')\n"
        "    audit_fix_dpmatter = audit_env('AUDIT_FIX_DPMATTER')\n"
        "    audit_fix_tensor = audit_env('AUDIT_FIX_TENSOR')\n"
        "    audit_polter_gr = audit_env('AUDIT_POLTER_GR')\n"
        "    audit_fix_metric_h = audit_env('AUDIT_FIX_METRIC_H')\n"
        "    audit_dgp_exact = audit_env('AUDIT_DGP_EXACT')\n"
        "    audit_fix_wde = audit_env('AUDIT_FIX_WDE')\n",
    ),
    (
        "    end subroutine GaugeInterface_Init\n",
        "    end subroutine GaugeInterface_Init\n\n"
        "    logical function audit_env(name)\n"
        "    character(len=*), intent(in) :: name\n"
        "    character(len=16) :: val\n"
        "    integer :: stat\n"
        "    call get_environment_variable(name, val, status=stat)\n"
        "    audit_env = (stat == 0 .and. len_trim(val) > 0 .and. val(1:1) /= '0')\n"
        "    end function audit_env\n",
    ),
    # 1: output double derivs call
    (
        "    call derivs(EV, EV%ScalEqsToPropagate, tau, y, yprime)\n\n    nullify(EV%OutputSources, EV%CustomSources)\n",
        "    call derivs(EV, EV%ScalEqsToPropagate, tau, y, yprime)\n"
        "    if (audit_fix_pidotdot) call derivs(EV, EV%ScalEqsToPropagate, tau, y, yprime)\n\n"
        "    nullify(EV%OutputSources, EV%CustomSources)\n",
    ),
    # 2: RSA dp_matter units
    (
        "            dp_matter = dgpnu\n",
        "            dp_matter = dgpnu\n            if (audit_fix_dpmatter) dp_matter = dgpnu/a2\n",
    ),
    # 3: tensor grhog_t/grhor_t
    (
        "    ! EFTCAMB MOD START: compute quantities for tensor equations\n",
        "    if (audit_fix_tensor) then\n        grhor_t = State%grhornomass/a2\n        grhog_t = State%grhog/a2\n"
        "    end if\n    ! EFTCAMB MOD START: compute quantities for tensor equations\n",
    ),
    # 4: original polterddot for EFTflag = 0
    (
        "            polterddot = pigddot/10._dl + 3._dl/5._dl*E2ddot\n            ! Original code:\n",
        "            polterddot = pigddot/10._dl + 3._dl/5._dl*E2ddot\n"
        "            if (audit_polter_gr .and. CP%EFTCAMB%EFTflag == 0) then\n"
        "                polterddot = -2.0d0/25.0d0*adotoa*dgq/(k*EV%Kf(1)) - 4.0d0/75.0d0*adotoa* &\n"
        "                    k*sigma - 4.0d0/75.0d0*dgpi - 2.0d0/75.0d0*dgrho/EV%Kf(1) - 3.0d0/ &\n"
        "                    50.0d0*k*octgdot*EV%Kf(2) + (1.0d0/25.0d0)*k*qgdot - 1.0d0/5.0d0 &\n"
        "                    *k*EV%Kf(2)*Edot(3) + (-1.0d0/10.0d0*pig + (7.0d0/10.0d0)* &\n"
        "                    polter - 3.0d0/5.0d0*E(2))*dopacity + (-1.0d0/10.0d0*pigdot &\n"
        "                    + (7.0d0/10.0d0)*polterdot - 3.0d0/5.0d0*Edot(2))*opacity\n"
        "            end if\n"
        "            ! Original code:\n",
    ),
    # 5: metric_h derivative while EFT inactive
    (
        "    else\n        ayprime(EV%w_ix) = 0._dl\n        ayprime(EV%w_ix + 1) = 0._dl\n    end if\n"
        "    ! EFTCAMB MOD END.\n\n    ! CDM equation of motion\n",
        "    else\n        ayprime(EV%w_ix) = 0._dl\n        ayprime(EV%w_ix + 1) = 0._dl\n"
        "        if (audit_fix_metric_h .and. CP%EFTCAMB%EFTCAMB_evolve_delta_phi .and. &\n"
        "            CP%EFTCAMB%EFTCAMB_evolve_metric_h) ayprime(EV%w_ix + 2) = 0._dl\n    end if\n"
        "    ! EFTCAMB MOD END.\n\n    ! CDM equation of motion\n",
    ),
    # 6: exact dgp
    (
        "    real(dl) grhormass_t, rhonu, qnu, clxnu, grhonu_t, gpnu_t, pnu\n",
        "    real(dl) grhormass_t, rhonu, qnu, clxnu, grhonu_t, gpnu_t, pnu\n"
        "    real(dl) audit_dpnu, audit_pinu, audit_clx, audit_q\n",
    ),
    (
        "            clxnu = y(EV%nu_ix(nu_i))\n            qnu = y(EV%nu_ix(nu_i) + 2)\n        else\n"
        "            ! Integrate over q\n            call Nu_Integrate_L01(EV, y, a, nu_i, clxnu, qnu)\n",
        "            clxnu = y(EV%nu_ix(nu_i))\n            qnu = y(EV%nu_ix(nu_i) + 2)\n"
        "            audit_dpnu = y(EV%nu_ix(nu_i) + 1)*rhonu\n        else\n"
        "            ! Integrate over q\n            call Nu_Integrate_L01(EV, y, a, nu_i, clxnu, qnu)\n"
        "            if (audit_dgp_exact) call Nu_Integrate_L012(EV, y, a, nu_i, audit_clx, audit_q, audit_dpnu, "
        "audit_pinu)\n",
    ),
    (
        "        if (present(dgp)) then\n            dgp = dgp + gpnu_t*clxnu\n        end if\n",
        "        if (present(dgp)) then\n            if (audit_dgp_exact) then\n"
        "                dgp = dgp + grhormass_t*audit_dpnu\n            else\n"
        "                dgp = dgp + gpnu_t*clxnu\n            end if\n        end if\n",
    ),
    # 7: w_dark_energy_t for EFT
    (
        "    else\n        grhov_t = 0._dl\n    end if\n    ! EFTCAMB MOD END.\n\n"
        "    ! EFTCAMB MOD START: start the total pressure",
        "    else\n        grhov_t = 0._dl\n        if (audit_fix_wde) w_dark_energy_t = -1._dl\n    end if\n"
        "    ! EFTCAMB MOD END.\n\n    ! EFTCAMB MOD START: start the total pressure",
    ),
]


def main():
    name = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else "release"
    src = ROOT / "trees" / "A0"
    dst = WORK / name
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, symlinks=True)
    eq = dst / "fortran" / "equations.f90"
    text = eq.read_text()
    for old, new in PATCHES:
        n = text.count(old)
        if n != 1:
            raise SystemExit(f"patch anchor found {n} times:\n{old}")
        text = text.replace(old, new)
    eq.write_text(text)
    (dst / "AUDIT_PATCHES.txt").write_text(__doc__)
    env = dict(os.environ, CAMB_BUILD_DIR=str(WORK / f"{name}_build"))
    log = WORK / f"build_{name}.log"
    if mode == "release":
        cmd = [PY, "setup.py", "make"]
        cwd = dst
    else:
        cmd = ["make", "python", "PYCAMB_OUTPUT_DIR=../camb/", f"FFLAGS={SNAN} {COMMON}"]
        cwd = dst / "fortran"
    with open(log, "w") as f:
        rc = subprocess.run(cmd, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT).returncode
    print(f"build {name} ({mode}) rc={rc} log={log}")


if __name__ == "__main__":
    main()

"""Build an instrumented copy of tree A0 for verifying finding evolve-metric-h-inconsistent (F-065).

Patches (both inactive unless the env var is set, so the build reproduces A0 otherwise):
  * equations.f90 derivs, delta_phi branch: if EFT_MH_LOG=<path>, append binary records (16 float64) for k within 3% of
    1e-3, 3e-3, 1e-2, 3e-2, 0.1, 1, 10 /Mpc:
      k, tau, a, h'(used), h'_00(constraint, recomputed here), h''(ii eq), etak, dgrho, dphi, dphi', RSA, TC, no_phot,
      dgq, etak', sigma
  * 09_EFTCAMB_IC.f90: if EFT_MH_ICFIX=1 and evolve_metric_h, start y(w_ix+2) from the EFT 00 constraint (delta_phi=0)
    instead of the GR h' = 2kz.

Usage: python metric_h_build_instr.py
"""

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

DEST = harness.AUDIT_ROOT / "work" / "verify" / "mh_instr"


def patch(path: Path, old: str, new: str) -> None:
    s = path.read_text()
    assert s.count(old) == 1, (path, old[:60], s.count(old))
    path.write_text(s.replace(old, new))


if not DEST.exists():
    shutil.copytree(harness.tree_path("A0"), DEST, symlinks=True)
    eq = DEST / "fortran" / "equations.f90"
    patch(
        eq,
        "    real(dl) ay(n), ayprime(n)\n    real(dl) tau, w\n    real(dl) k, k2\n",
        "    real(dl) ay(n), ayprime(n)\n    real(dl) tau, w\n    real(dl) k, k2\n"
        "    integer, save :: mhl_unit = -2, mhl_count = 0\n"
        "    character(len=512) :: mhl_path\n"
        "    integer :: mhl_st\n"
        "    real(dl) :: mhl_hc\n"
        "    real(dl), parameter :: mhl_k(7) = [1.e-3_dl, 3.e-3_dl, 1.e-2_dl, 3.e-2_dl, 0.1_dl, 1._dl, 10._dl]\n",
    )
    patch(
        eq,
        "            sigma = (metric_h_dot + 6._dl*ayprime(ix_etak)/k)/(2._dl*k)\n"
        "            EV%eft_cache%etak = etak\n"
        "            EV%eft_cache%etakdot = ayprime(ix_etak)\n"
        "            EV%eft_cache%sigma = sigma\n",
        "            sigma = (metric_h_dot + 6._dl*ayprime(ix_etak)/k)/(2._dl*k)\n"
        "            EV%eft_cache%etak = etak\n"
        "            EV%eft_cache%etakdot = ayprime(ix_etak)\n"
        "            EV%eft_cache%sigma = sigma\n"
        "            if (mhl_unit == -2) then\n"
        "                !$OMP CRITICAL (mhlog_init)\n"
        "                if (mhl_unit == -2) then\n"
        "                    call get_environment_variable('EFT_MH_LOG', mhl_path, status=mhl_st)\n"
        "                    if (mhl_st == 0 .and. len_trim(mhl_path) > 0) then\n"
        "                        open(newunit=mhl_unit, file=trim(mhl_path), access='stream', form='unformatted', &\n"
        "                            status='replace')\n"
        "                    else\n"
        "                        mhl_unit = -1\n"
        "                    end if\n"
        "                end if\n"
        "                !$OMP END CRITICAL (mhlog_init)\n"
        "            end if\n"
        "            if (mhl_unit /= -1 .and. mhl_unit /= -2) then\n"
        "                if (any(abs(k/mhl_k - 1._dl) < 0.03_dl)) then\n"
        "                    mhl_hc = (4._dl*k*etak + 2._dl*dgrho/Meff2)/(2._dl - alphaB)/adotoa + &\n"
        "                        s00p*delta_phi_dot + (s00 + s00k*k2)*delta_phi\n"
        "                    !$OMP CRITICAL (mhlog)\n"
        "                    write(mhl_unit) k, tau, a, metric_h_dot, mhl_hc, metric_h_ddot, etak, dgrho, delta_phi, &\n"
        "                        delta_phi_dot, merge(1._dl, 0._dl, EV%no_nu_multpoles), &\n"
        "                        merge(1._dl, 0._dl, EV%TightCoupling), merge(1._dl, 0._dl, EV%no_phot_multpoles), &\n"
        "                        dgq, ayprime(ix_etak), sigma\n"
        "                    mhl_count = mhl_count + 1\n"
        "                    if (mod(mhl_count, 50) == 0) flush(mhl_unit)\n"
        "                    !$OMP END CRITICAL (mhlog)\n"
        "                end if\n"
        "            end if\n",
    )
    ic = DEST / "fortran" / "eftcamb" / "09_EFTCAMB_IC.f90"
    patch(
        ic,
        "    real(dl) :: wnu_arr(max_nu)\n",
        "    real(dl) :: wnu_arr(max_nu)\n    character(len=16) :: mhl_env\n    integer :: mhl_st\n",
    )
    patch(
        ic,
        "        if (CP%EFTCAMB%EFTCAMB_evolve_metric_h) y(EV%w_ix + 2) = 2._dl*k*z\n",
        "        if (CP%EFTCAMB%EFTCAMB_evolve_metric_h) then\n"
        "            call get_environment_variable('EFT_MH_ICFIX', mhl_env, status=mhl_st)\n"
        "            if (mhl_st == 0 .and. trim(mhl_env) == '1') then\n"
        "                y(EV%w_ix + 2) = (4._dl*k*etak + 2._dl*dgrho/EV%eft_cache%Meff2)/ &\n"
        "                    (2._dl + 2._dl*EV%eft_cache%alphaB)/adotoa\n"
        "            else\n"
        "                y(EV%w_ix + 2) = 2._dl*k*z\n"
        "            end if\n"
        "        end if\n",
    )
    print("patched", DEST)
print(harness.build_tree(DEST, DEST))

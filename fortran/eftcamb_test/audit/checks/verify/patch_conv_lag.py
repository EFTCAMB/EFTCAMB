"""Verifier instrumentation for finding fullmap-bg-eft-stale-background (independent of t2d-static's LAG patch).

Copies tree A0 to work/verify/lag/CONV and patches fortran/equations.f90 (derivs, full models only): after the
second compute_H_derivs call, it iterates the whole full-model background block (bg EFT functions -> compute_adotoa ->
massive-nu first derivatives -> compute_H_derivs -> massive-nu second derivatives -> compute_H_derivs) on a COPY of
the cache until it reaches the self-consistent fixed point (all of adotoa, Hdot, EFTc, EFTLambda stationary to 1e-12),
and compares the values derivs actually used with the fixed point, normalised by natural scales:
  1 |H_used/H_conv - 1|          2 |dHdot|/H^2          3 |dEFTc|/(3H^2)          4 |dEFTLambda|/(3H^2)
  5 |dEFTcdot|/(3H^3)           6 |dEFTLambdadot|/(3H^3) 7 |EFTc_used/EFTc_conv-1| 8 |dHdotdot|/H^3
Running maxima (a >= 1e-3) are written to stderr as CONVMAX lines. With EFT_CONV_FIX=1 derivs continues with the
fixed-point cache (full fix of the lag, including the BH self-referential adotoa); otherwise the used cache is kept.
Single-thread instrumentation (module state); run with OMP_NUM_THREADS=1.

Usage: python patch_conv_lag.py   (then build with harness.build_tree(A0, dest))
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
SRC = AUDIT / "trees" / "A0"
DEST = AUDIT / "work" / "verify" / "lag" / "CONV"

MODVARS = """    integer, parameter :: basic_num_eqns = 4
    ! AUDIT CONV instrumentation (module state, single thread only)
    real(dl), save :: conv_max(8) = 0._dl
    integer, save :: conv_calls = 0, conv_itmax = 0, conv_nonconv = 0
    logical, save :: conv_init = .false., conv_fix = .false.
"""

DECL = """    real(dl) xe, Trad, Delta_TM, Tmat, Delta_TCMB
    type(TEFTCAMB_timestep_cache) :: conv_c
    real(dl) conv_prev(4), conv_new(4), conv_used(8), conv_fp(8), conv_m(8), conv_H, conv_sc(4)
    integer conv_it, conv_i
    character(len=8) conv_env
"""

NU1 = """                conv_c%grhonu_tot = 0._dl
                conv_c%gpinu_tot = 0._dl
                conv_c%grhonudot_tot = 0._dl
                conv_c%gpinudot_tot = 0._dl
                if (CP%Num_Nu_massive /= 0) then
                    do nu_i = 1, CP%Nu_mass_eigenstates
                        EFT_grhonu = 0._dl
                        EFT_gpinu = 0._dl
                        grhormass_t = State%grhormass(nu_i)/a**2
                        call ThermalNuBackground%rho_P(a*State%nu_masses(nu_i), EFT_grhonu, EFT_gpinu)
                        conv_c%grhonu_tot = conv_c%grhonu_tot + grhormass_t*EFT_grhonu
                        conv_c%gpinu_tot = conv_c%gpinu_tot + grhormass_t*EFT_gpinu
                        conv_c%gpinudot_tot = conv_c%gpinudot_tot + &
                            grhormass_t*(ThermalNuBackground%pidot(a*State%nu_masses(nu_i), conv_H, EFT_gpinu) - &
                            4._dl*conv_H*EFT_gpinu)
                        conv_c%grhonudot_tot = conv_c%grhonudot_tot + &
                            grhormass_t*(ThermalNuBackground%drho(a*State%nu_masses(nu_i), conv_H) - &
                            4._dl*conv_H*EFT_grhonu)
                    end do
                end if
                conv_c%gpresdotm_t = -4._dl*conv_H*(grhog_t + grhor_t)/3._dl + conv_c%gpinudot_tot
                call CP%EFTCAMB%model%compute_H_derivs(a, CP%eft_par_cache, conv_c)
"""

NU2 = """                conv_c%gpinudotdot_tot = 0._dl
                if (CP%Num_Nu_massive /= 0) then
                    do nu_i = 1, CP%Nu_mass_eigenstates
                        EFT_grhonu = 0._dl
                        EFT_gpinu = 0._dl
                        grhormass_t = State%grhormass(nu_i)/a**2
                        call ThermalNuBackground%rho_P(a*State%nu_masses(nu_i), EFT_grhonu, EFT_gpinu)
                        EFT_gpinudot = ThermalNuBackground%pidot(a*State%nu_masses(nu_i), conv_H, EFT_gpinu)
                        conv_c%gpinudotdot_tot = conv_c%gpinudotdot_tot - &
                            4._dl*conv_H*grhormass_t*(EFT_gpinudot - 4._dl*conv_H*EFT_gpinu) + &
                            grhormass_t*(ThermalNuBackground%pidotdot(a*State%nu_masses(nu_i), conv_H, conv_c%Hdot, &
                            EFT_gpinu, EFT_gpinudot) - 4._dl*conv_c%Hdot*EFT_gpinu - 4._dl*conv_H*EFT_gpinudot)
                    end do
                end if
                conv_c%gpresdotdotm_t = -(4._dl/3._dl)*conv_c%Hdot*(grhog_t + grhor_t) + &
                    (16._dl/3._dl)*(grhog_t + grhor_t)*conv_H**2 + conv_c%gpinudotdot_tot
                call CP%EFTCAMB%model%compute_H_derivs(a, CP%eft_par_cache, conv_c)
"""

BLOCK = """        ! compute third H derivate:
        call CP%EFTCAMB%model%compute_H_derivs(a, CP%eft_par_cache, EV%eft_cache)
        ! AUDIT CONV instrumentation:
        if (.not. CP%EFTCAMB%EFTCAMB_model_is_designer) then
            if (.not. conv_init) then
                conv_env = ''
                call get_environment_variable('EFT_CONV_FIX', conv_env)
                conv_fix = trim(conv_env) == '1'
                conv_init = .true.
                write(0, '(a,l2)') 'CONVINIT fix =', conv_fix
            end if
            conv_c = EV%eft_cache
            do conv_it = 1, 400
                conv_prev = [conv_c%adotoa, conv_c%Hdot, conv_c%EFTc, conv_c%EFTLambda]
                call CP%EFTCAMB%model%compute_background_EFT_functions(a, CP%eft_par_cache, conv_c)
                call CP%EFTCAMB%model%compute_adotoa(a, CP%eft_par_cache, conv_c)
                conv_H = conv_c%adotoa
""" + NU1 + NU2 + """                conv_new = [conv_c%adotoa, conv_c%Hdot, conv_c%EFTc, conv_c%EFTLambda]
                conv_sc = [conv_c%adotoa, conv_c%adotoa**2, 3._dl*conv_c%adotoa**2, 3._dl*conv_c%adotoa**2]
                if (all(abs(conv_new - conv_prev) <= 1.d-12*conv_sc)) exit
            end do
            conv_calls = conv_calls + 1
            if (conv_it > 400) conv_nonconv = conv_nonconv + 1
            conv_itmax = max(conv_itmax, min(conv_it, 400))
            conv_H = conv_c%adotoa
            conv_used = [EV%eft_cache%adotoa, EV%eft_cache%Hdot, EV%eft_cache%EFTc, EV%eft_cache%EFTLambda, &
                EV%eft_cache%EFTcdot, EV%eft_cache%EFTLambdadot, EV%eft_cache%EFTc, EV%eft_cache%Hdotdot]
            conv_fp = [conv_c%adotoa, conv_c%Hdot, conv_c%EFTc, conv_c%EFTLambda, &
                conv_c%EFTcdot, conv_c%EFTLambdadot, conv_c%EFTc, conv_c%Hdotdot]
            conv_m(1) = abs(conv_used(1)/conv_fp(1) - 1._dl)
            conv_m(2) = abs(conv_used(2) - conv_fp(2))/conv_H**2
            conv_m(3) = abs(conv_used(3) - conv_fp(3))/(3._dl*conv_H**2)
            conv_m(4) = abs(conv_used(4) - conv_fp(4))/(3._dl*conv_H**2)
            conv_m(5) = abs(conv_used(5) - conv_fp(5))/(3._dl*conv_H**3)
            conv_m(6) = abs(conv_used(6) - conv_fp(6))/(3._dl*conv_H**3)
            conv_m(7) = abs(conv_used(7) - conv_fp(7))/max(abs(conv_fp(7)), 1.d-300)
            conv_m(8) = abs(conv_used(8) - conv_fp(8))/conv_H**3
            if (a >= 1.d-3) then
                do conv_i = 1, 8
                    if (conv_m(conv_i) > conv_max(conv_i)) then
                        conv_max(conv_i) = conv_m(conv_i)
                        write(0, '(a,i2,i4,7es14.5)') 'CONVMAX ', conv_i, min(conv_it, 400), a, k, tau, &
                            conv_used(conv_i), conv_fp(conv_i), conv_m(conv_i), conv_H
                    end if
                end do
            end if
            if (mod(conv_calls, 200000) == 0) write(0, '(a,i10,2i5,8es11.3)') 'CONVSUM ', conv_calls, conv_itmax, &
                conv_nonconv, conv_max
            if (conv_fix) then
                EV%eft_cache = conv_c
                adotoa = conv_c%adotoa
            end if
        end if
"""

ANCHOR_BLOCK = """        ! compute third H derivate:
        call CP%EFTCAMB%model%compute_H_derivs(a, CP%eft_par_cache, EV%eft_cache)
"""


def patch(text: str) -> str:
    for old, new in [("    integer, parameter :: basic_num_eqns = 4\n", MODVARS),
                     ("    real(dl) xe, Trad, Delta_TM, Tmat, Delta_TCMB\n", DECL),
                     (ANCHOR_BLOCK, BLOCK)]:
        n = text.count(old)
        if n != 1:
            sys.exit(f"anchor found {n} times: {old!r}")
        text = text.replace(old, new)
    return text


def main() -> None:
    if not DEST.exists():
        shutil.copytree(SRC, DEST, symlinks=True)
    eq = DEST / "fortran" / "equations.f90"
    orig = (SRC / "fortran" / "equations.f90").read_text()
    eq.write_text(patch(orig))
    print(f"patched {eq}")


if __name__ == "__main__":
    main()

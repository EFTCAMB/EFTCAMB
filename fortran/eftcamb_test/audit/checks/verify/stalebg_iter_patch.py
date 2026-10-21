"""Verifier (fullmap-bg-eft-stale-background): independent fix-by-iteration tree.

After the second compute_H_derivs in derivs, for full models, iterate (bg EFT functions -> adotoa -> H derivs)
EFT_ITER times (env var, default 0 = code under audit unchanged), then refresh the local adotoa. Prints to stderr the
max relative change of EFTc / EFTLambda / adotoa between the value derivs had and the iterated value (a >= 1e-3).
Usage: python stalebg_iter_patch.py <tree copy>   (patches <tree>/fortran/equations.f90 in place)
"""

import sys
from pathlib import Path

f = Path(sys.argv[1]) / "fortran" / "equations.f90"
s = f.read_text()
anchor = """        call CP%EFTCAMB%model%compute_H_derivs(a, CP%eft_par_cache, EV%eft_cache)
        ! store:
        adotdota = EV%eft_cache%Hdot + EV%eft_cache%adotoa**2"""
assert s.count(anchor) == 1
new = """        call CP%EFTCAMB%model%compute_H_derivs(a, CP%eft_par_cache, EV%eft_cache)
        ! VERIFY ITER instrumentation
        if (.not. CP%EFTCAMB%EFTCAMB_model_is_designer) then
            if (.not. vit_init) then
                vit_env = ''
                call get_environment_variable('EFT_ITER', vit_env)
                if (len_trim(vit_env) > 0) read(vit_env, *) vit_n
                vit_init = .true.
                write(0, '(a,i4)') 'VITINIT n =', vit_n
            end if
            if (vit_n > 0) then
                vit_old = [EV%eft_cache%EFTc, EV%eft_cache%EFTLambda, EV%eft_cache%adotoa]
                do vit_i = 1, vit_n
                    call CP%EFTCAMB%model%compute_background_EFT_functions(a, CP%eft_par_cache, EV%eft_cache)
                    call CP%EFTCAMB%model%compute_adotoa(a, CP%eft_par_cache, EV%eft_cache)
                    call CP%EFTCAMB%model%compute_H_derivs(a, CP%eft_par_cache, EV%eft_cache)
                end do
                call CP%EFTCAMB%model%compute_background_EFT_functions(a, CP%eft_par_cache, EV%eft_cache)
                adotoa = EV%eft_cache%adotoa
                vit_new = [EV%eft_cache%EFTc, EV%eft_cache%EFTLambda, EV%eft_cache%adotoa]
                do vit_i = 1, 3
                    vit_rel = abs(vit_old(vit_i) - vit_new(vit_i))/max(abs(vit_new(vit_i)), 1.d-300)
                    if (a >= 1.d-3 .and. vit_rel > vit_max(vit_i)) then
                        vit_max(vit_i) = vit_rel
                        write(0, '(a,i2,5es14.5)') 'VITMAX ', vit_i, a, k, vit_old(vit_i), vit_new(vit_i), vit_rel
                    end if
                end do
            end if
        end if
        ! store:
        adotdota = EV%eft_cache%Hdot + EV%eft_cache%adotoa**2"""
s = s.replace(anchor, new)
a2 = "    integer, parameter :: basic_num_eqns = 4\n"
assert s.count(a2) == 1
s = s.replace(a2, a2 + "    real(dl), save :: vit_max(3) = 0._dl\n    integer, save :: vit_n = 0\n"
              "    logical, save :: vit_init = .false.\n")
a3 = "    real(dl) xe, Trad, Delta_TM, Tmat, Delta_TCMB\n"
assert s.count(a3) == 1
s = s.replace(a3, a3 + "    real(dl) vit_old(3), vit_new(3), vit_rel\n    integer vit_i\n    character(len=16) vit_env\n", 1)
f.write_text(s)
print("patched", f)

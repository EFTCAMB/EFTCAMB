"""Patch a copy of A0 with environment-variable switches for attributing the EFT_IC_type 2/3 problem.

IC23_PI_MODE   = series (default, unchanged) | attractor (use the EFT_IC_type=1 pi, pi' formula) | zero
IC23_ADIAB     = mg (default, unchanged) | gr (use the GR adiabatic CAMB ICs even when EFT_IC_type /= 1)
IC23_DUMP      = 1 -> print k, tau, a, k*tau, pi_series, pi_attractor, pi'_series, pi'_attractor at turn-on (stderr)
"""

import sys
from pathlib import Path

tree = Path(sys.argv[1])
ic = tree / "fortran/eftcamb/09_EFTCAMB_IC.f90"
s = ic.read_text()

old_decl = "    real(dl) :: Omega0, Gamma10, Gamma20, Gamma30\n    real(dl) :: wnu_arr(max_nu)\n"
assert s.count(old_decl) == 1
s = s.replace(old_decl, old_decl + "    character(len=32) :: sw_mode, sw_dump\n    integer :: sw_st, sw_st2\n"
              "    real(dl) :: sw_pa, sw_dpa, sw_den\n")

old_end = """        case default
            y(EV%w_ix) = 0._dl
            y(EV%w_ix + 1) = 0._dl
        end select
    end if

    ! reset the cache:"""
assert s.count(old_end) == 1
new_end = """        case default
            y(EV%w_ix) = 0._dl
            y(EV%w_ix + 1) = 0._dl
        end select
        ! --- verify instrumentation ---
        sw_den = EV%eft_cache%EFTpiC + k2*EV%eft_cache%EFTpiD1 + k2**2*EV%eft_cache%EFTpiD2
        if (sw_den /= 0._dl) then
            sw_pa = -CP%eft_par_cache%h0_Mpc*EV%eft_cache%EFTpiE/sw_den
            sw_dpa = CP%eft_par_cache%h0_Mpc*(EV%eft_cache%EFTpiE/sw_den**2*a*EV%eft_cache%adotoa* &
                (EFTpiCdotFunction(a, k) + k*k*EFTpiDdotFunction(a, k)) - EFTpiEdot/sw_den)
        else
            sw_pa = 0._dl
            sw_dpa = 0._dl
        end if
        call get_environment_variable('IC23_DUMP', sw_dump, status=sw_st2)
        if (sw_st2 == 0) then
            write(0, '(a,i2,9es15.6)') 'IC23DUMP ', CP%EFTCAMB%EFT_IC_type, k, tau, a, k*tau, y(EV%w_ix), sw_pa, &
                y(EV%w_ix + 1), sw_dpa, EV%eft_cache%EFTOmegaV
            flush(0)
        end if
        call get_environment_variable('IC23_PI_MODE', sw_mode, status=sw_st)
        if (sw_st == 0 .and. CP%EFTCAMB%EFT_IC_type /= 1) then
            if (trim(sw_mode) == 'attractor') then
                y(EV%w_ix) = sw_pa
                y(EV%w_ix + 1) = sw_dpa
            else if (trim(sw_mode) == 'zero') then
                y(EV%w_ix) = 0._dl
                y(EV%w_ix + 1) = 0._dl
            end if
        end if
    end if

    ! reset the cache:"""
s = s.replace(old_end, new_end)
ic.write_text(s)

eq = tree / "fortran/equations.f90"
t = eq.read_text()
old = "    real(dl) :: a_eval\n"
assert t.count(old) == 1
t = t.replace(old, old + "    character(len=32) :: sw_adiab\n    integer :: sw_sta\n    logical :: sw_gr\n")
old2 = """    chi = 1 ! Get transfer function for chi
    ! EFTCAMB MOD START
    if (CP%EFTCAMB%EFTflag == 0 .or. CP%EFTCAMB%EFT_IC_type == 1) then"""
assert t.count(old2) == 1
t = t.replace(old2, """    chi = 1 ! Get transfer function for chi
    ! EFTCAMB MOD START
    call get_environment_variable('IC23_ADIAB', sw_adiab, status=sw_sta)
    sw_gr = (sw_sta == 0 .and. trim(sw_adiab) == 'gr')
    if (CP%EFTCAMB%EFTflag == 0 .or. CP%EFTCAMB%EFT_IC_type == 1 .or. sw_gr) then""")
eq.write_text(t)
print("patched", tree)

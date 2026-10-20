"""Source patches applied to audit tree COPIES (never to the repo). Each patch is an exact old->new replacement.

Usage: python patch_tools.py TREE PATCH [PATCH ...]
Patches:
  f002            upstream results.f90 size() of unallocated array (nested if)
  derivst_grho    equations.f90 derivst: use EV%eft_cache%grhog_t/grhor_t (locals are used before assignment)
  bh_hpprime      008p5 Beyond Horndeski: drop the dead x3ppp statement that reads uninitialised hpprime
  audit_mass      09_EFTCAMB_stability: print analytic vs numerical additional derivatives and mass eigenvalues
  fix_desaddder   06p2 designer compute_additional_derivs: evaluate helper at temp_a with temp_cache
  fix_fulladdder  06p1 full compute_additional_derivs: evaluate helper at temp_a with temp_cache, set all outputs
  w_de_init, ev_reset, printing_reset, printing_probe, quint_drag_y2, full_bg_recompute, interp_clamp,
  cli_p_eftcamb, qs_printing_only: see the comments above each PATCHES entry
"""

from __future__ import annotations

import sys
from pathlib import Path

PATCHES: dict[str, list[tuple[str, str, str]]] = {}

PATCHES["f002"] = [(
    "fortran/results.f90",
    """    if (allocated(this%tb) .and. this%nthermo /= size(this%tb)) then
        deallocate(this%ScaleFactor, this%cs2, this%ddotmu)
        deallocate(this%adot, this%dadot)
        deallocate(this%tb, this%xe, this%emmu, this%dotmu)
        deallocate(this%demmu, this%dddotmu, this%ddddotmu)
        deallocate(this%thermo_horner)
        if (dowinlens .and. allocated(this%winlens)) deallocate(this%winlens, this%dwinlens)
    end if""",
    """    if (allocated(this%tb)) then
    if (this%nthermo /= size(this%tb)) then
        deallocate(this%ScaleFactor, this%cs2, this%ddotmu)
        deallocate(this%adot, this%dadot)
        deallocate(this%tb, this%xe, this%emmu, this%dotmu)
        deallocate(this%demmu, this%dddotmu, this%ddddotmu)
        deallocate(this%thermo_horner)
        if (dowinlens .and. allocated(this%winlens)) deallocate(this%winlens, this%dwinlens)
    end if
    end if""")]

PATCHES["derivst_grho"] = [(
    "fortran/equations.f90",
    """            EV%eft_cache%gpresdotdotm_t = -(4._dl/3._dl)*EV%eft_cache%Hdot*(grhog_t + grhor_t) + &
                (16._dl/3._dl)*(grhog_t + grhor_t)*adotoa**2 + EV%eft_cache%gpinudotdot_tot""",
    """            EV%eft_cache%gpresdotdotm_t = -(4._dl/3._dl)*EV%eft_cache%Hdot*(EV%eft_cache%grhog_t + EV%eft_cache%grhor_t) + &
                (16._dl/3._dl)*(EV%eft_cache%grhog_t + EV%eft_cache%grhor_t)*adotoa**2 + EV%eft_cache%gpinudotdot_tot""")]

PATCHES["bh_hpprime"] = [(
    "fortran/eftcamb/08f_full_models/008p5_Beyond_Horndeski.f90",
    """    x4ppp = 4._dl*(x4pp*epsilon + 2._dl*x4p*epsilon_p + x4*epsilon_pp)
    x3ppp = -x3*(hpprime - 3._dl*epsilon_pp) - 2._dl*x3p*(hprime - 3._dl*epsilon_p) - (h - 3._dl*epsilon)*x3pp
    call self%compute_adotoa(a, eft_par_cache, eft_cache)""",
    """    x4ppp = 4._dl*(x4pp*epsilon + 2._dl*x4p*epsilon_p + x4*epsilon_pp)
    call self%compute_adotoa(a, eft_par_cache, eft_cache)""")]

_MASS_OLD = """    if (a > input_EFTCAMB%EFTCAMB_pert_turn_on .and. input_EFTCAMB%EFT_mass_stability) then
        call input_model%compute_additional_derivs(a, params_cache, eft_cache)
        call input_model%compute_mass_factors(a, params_cache, eft_cache)
    end if"""
_MASS_NEW = """    if (a > input_EFTCAMB%EFTCAMB_pert_turn_on .and. input_EFTCAMB%EFT_mass_stability) then
        audit_cache = eft_cache
        call input_model%compute_mass_factors(a, params_cache, audit_cache)
        write(*, '(a,12ES17.8)') 'AUDIT_ANALYTIC ', a, eft_cache%EFTcdotdot, eft_cache%EFTLambdadotdot, &
            eft_cache%EFTGamma1PP, eft_cache%EFTGamma2PP, eft_cache%EFTGamma3PP, eft_cache%EFTGamma2PPP, &
            eft_cache%EFTGamma3PPP, eft_cache%EFTGamma3PPPP, audit_cache%EFT_mu1, audit_cache%EFT_mu2, eft_cache%adotoa
        call input_model%compute_additional_derivs(a, params_cache, eft_cache)
        call input_model%compute_mass_factors(a, params_cache, eft_cache)
        write(*, '(a,12ES17.8)') 'AUDIT_NUMERIC  ', a, eft_cache%EFTcdotdot, eft_cache%EFTLambdadotdot, &
            eft_cache%EFTGamma1PP, eft_cache%EFTGamma2PP, eft_cache%EFTGamma3PP, eft_cache%EFTGamma2PPP, &
            eft_cache%EFTGamma3PPP, eft_cache%EFTGamma3PPPP, eft_cache%EFT_mu1, eft_cache%EFT_mu2, eft_cache%adotoa
    end if"""
PATCHES["audit_mass"] = [
    ("fortran/eftcamb/09_EFTCAMB_stability.f90", _MASS_OLD, _MASS_NEW),
    ("fortran/eftcamb/09_EFTCAMB_stability.f90",
     """    real(dl) :: grhonu, gpinu, grhormass_t, grhonudot, gpinudot, gpinudotdot
    integer  :: nu_i, ind, ind_max
""",
     """    real(dl) :: grhonu, gpinu, grhormass_t, grhonudot, gpinudot, gpinudotdot
    integer  :: nu_i, ind, ind_max
    type(TEFTCAMB_timestep_cache) :: audit_cache
""")]

_DES = "fortran/eftcamb/06p2_abstract_EFTCAMB_designer.f90"
PATCHES["fix_desaddder"] = [
    (_DES, """    call self%compute_adotoa(a, temp_par_cache, eft_cache)
    ! store adotoa:
    adotoa = temp_cache%adotoa""", """    call self%compute_adotoa(temp_a, temp_par_cache, temp_cache)
    ! store adotoa:
    adotoa = temp_cache%adotoa"""),
    (_DES, """    cP = eft_cache%EFTcdot/(a*eft_cache%adotoa) + 2._dl*eft_cache%EFTc/a

    end subroutine dummy_helper""", """    cP = temp_cache%EFTcdot/(temp_a*temp_cache%adotoa) + 2._dl*temp_cache%EFTc/temp_a

    end subroutine dummy_helper"""),
]
# the neutrino loops of the designer helper use the host a and host eft_cache%Hdot
for _old, _new in (
    ("            grhormass_t = temp_par_cache%grhormass(nu_i)/a**2\n"
     "            call ThermalNuBackground%rho_P(a*temp_par_cache%nu_masses(nu_i), EFT_grhonu, EFT_gpinu)",
     "            grhormass_t = temp_par_cache%grhormass(nu_i)/temp_a**2\n"
     "            call ThermalNuBackground%rho_P(temp_a*temp_par_cache%nu_masses(nu_i), EFT_grhonu, EFT_gpinu)"),
):
    PATCHES["fix_desaddder"].append((_DES, _old, _new, "all"))  # type: ignore[arg-type]
for _old, _new in (
    ("ThermalNuBackground%drho(a*temp_par_cache%nu_masses(nu_i), adotoa)",
     "ThermalNuBackground%drho(temp_a*temp_par_cache%nu_masses(nu_i), adotoa)"),
    ("ThermalNuBackground%pidot(a*temp_par_cache%nu_masses(nu_i), adotoa, EFT_gpinu)",
     "ThermalNuBackground%pidot(temp_a*temp_par_cache%nu_masses(nu_i), adotoa, EFT_gpinu)"),
    ("ThermalNuBackground%pidotdot(a*temp_par_cache%nu_masses(nu_i), adotoa, eft_cache%Hdot, &",
     "ThermalNuBackground%pidotdot(temp_a*temp_par_cache%nu_masses(nu_i), adotoa, temp_cache%Hdot, &"),
    ("EFT_gpinu, EFT_gpinudot) - 4._dl*eft_cache%Hdot*EFT_gpinu - 4._dl*adotoa*EFT_gpinudot)",
     "EFT_gpinu, EFT_gpinudot) - 4._dl*temp_cache%Hdot*EFT_gpinu - 4._dl*adotoa*EFT_gpinudot)"),
):
    PATCHES["fix_desaddder"].append((_DES, _old, _new, "all"))  # type: ignore[arg-type]

_FULL = "fortran/eftcamb/06p1_abstract_EFTCAMB_full_map.f90"
PATCHES["fix_fulladdder"] = [
    (_FULL, """    call self%compute_background_EFT_functions(a, temp_par_cache, eft_cache)
    call self%compute_adotoa(a, temp_par_cache, eft_cache)""",
     """    call self%compute_background_EFT_functions(temp_a, temp_par_cache, temp_cache)
    call self%compute_adotoa(temp_a, temp_par_cache, temp_cache)"""),
    (_FULL, """    gamma3P = temp_cache%EFTGamma3P
    cP = eft_cache%EFTcdot/(a*eft_cache%adotoa) + 2._dl*eft_cache%EFTc/a
    LambdaP = eft_cache%EFTLambdadot/(a*eft_cache%adotoa) + 2._dl*eft_cache%EFTLambda/a""",
     """    gamma3P = temp_cache%EFTGamma3P
    gamma2PP = temp_cache%EFTGamma2PP
    gamma3PP = temp_cache%EFTGamma3PP
    gamma3PPP = temp_cache%EFTGamma3PPP
    cP = temp_cache%EFTcdot/(temp_a*temp_cache%adotoa) + 2._dl*temp_cache%EFTc/temp_a
    LambdaP = temp_cache%EFTLambdadot/(temp_a*temp_cache%adotoa) + 2._dl*temp_cache%EFTLambda/temp_a"""),
]
for _old, _new in (
    ("            grhormass_t = temp_par_cache%grhormass(nu_i)/a**2\n"
     "            call ThermalNuBackground%rho_P(a*temp_par_cache%nu_masses(nu_i), EFT_grhonu, EFT_gpinu)",
     "            grhormass_t = temp_par_cache%grhormass(nu_i)/temp_a**2\n"
     "            call ThermalNuBackground%rho_P(temp_a*temp_par_cache%nu_masses(nu_i), EFT_grhonu, EFT_gpinu)"),
    ("ThermalNuBackground%drho(a*temp_par_cache%nu_masses(nu_i), adotoa)",
     "ThermalNuBackground%drho(temp_a*temp_par_cache%nu_masses(nu_i), adotoa)"),
    ("ThermalNuBackground%pidot(a*temp_par_cache%nu_masses(nu_i), adotoa, EFT_gpinu)",
     "ThermalNuBackground%pidot(temp_a*temp_par_cache%nu_masses(nu_i), adotoa, EFT_gpinu)"),
    ("ThermalNuBackground%pidotdot(a*temp_par_cache%nu_masses(nu_i), adotoa, eft_cache%Hdot, &",
     "ThermalNuBackground%pidotdot(temp_a*temp_par_cache%nu_masses(nu_i), adotoa, temp_cache%Hdot, &"),
    ("EFT_gpinu, EFT_gpinudot) - 4._dl*eft_cache%Hdot*EFT_gpinu - 4._dl*adotoa*EFT_gpinudot)",
     "EFT_gpinu, EFT_gpinudot) - 4._dl*temp_cache%Hdot*EFT_gpinu - 4._dl*adotoa*EFT_gpinudot)"),
):
    PATCHES["fix_fulladdder"].append((_FULL, _old, _new, "all"))  # type: ignore[arg-type]


PATCHES["w_de_init"] = [(
    "fortran/equations.f90",
    """    else
        grhov_t = 0._dl
    end if
    ! EFTCAMB MOD END.""",
    """    else
        grhov_t = 0._dl
        w_dark_energy_t = -1._dl
    end if
    ! EFTCAMB MOD END.""")]

PATCHES["ev_reset"] = [(
    "fortran/cmbmain.f90",
    """    EV%q = ThisSources%Evolve_q%points(q_ix)
""",
    """    EV%q = ThisSources%Evolve_q%points(q_ix)
    EV%EFTCAMB_printing = .false.
    EV%EFTCAMB_QS = .false.
    call EV%eft_cache%initialize()
""")]

PATCHES["printing_reset"] = [(
    "fortran/cmbmain.f90",
    """    EV%q = ThisSources%Evolve_q%points(q_ix)
""",
    """    EV%q = ThisSources%Evolve_q%points(q_ix)
    EV%EFTCAMB_printing = .false.
    EV%EFTCAMB_QS = .false.
""")]

# instrumentation: report (once per k mode) when the OMP-private EV arrives with EFTCAMB_printing = .true.
PATCHES["printing_probe"] = [(
    "fortran/cmbmain.f90",
    """    EV%q = ThisSources%Evolve_q%points(q_ix)
""",
    """    EV%q = ThisSources%Evolve_q%points(q_ix)
    if (EV%EFTCAMB_printing) write(*, '(a,i6)') 'AUDIT_PRINTING_TRUE q_ix=', q_ix
    if (.not. EV%EFTCAMB_printing) write(*, '(a,i6)') 'AUDIT_PRINTING_FALSE q_ix=', q_ix
""")]

PATCHES["quint_drag_y2"] = [(
    "fortran/eftcamb/08f_full_models/008p4_quintessence.f90",
    """    y(1) = phi_ini
""",
    """    y(1) = phi_ini
    y(2) = 0._dl
""")]

# recompute the background EFT functions of full-mapping models after H and its derivatives are known at this time
_RECOMP_OLD = """        ! finish computing the background:
        if (CP%EFTCAMB%EFTCAMB_model_is_designer) then
            call CP%EFTCAMB%model%compute_background_EFT_functions(a, CP%eft_par_cache, EV%eft_cache)
        end if"""
_RECOMP_NEW = """        ! finish computing the background:
        call CP%EFTCAMB%model%compute_background_EFT_functions(a, CP%eft_par_cache, EV%eft_cache)"""
PATCHES["full_bg_recompute"] = [
    ("fortran/equations.f90", _RECOMP_OLD, _RECOMP_NEW),
    ("fortran/equations.f90", "\n".join("    " + ln for ln in _RECOMP_OLD.splitlines()),
     "\n".join("    " + ln for ln in _RECOMP_NEW.splitlines())),
]


# clamp the equispaced-interpolation index so that x within rounding of x_final cannot address x(N+1), y(N+1)
PATCHES["interp_clamp"] = [(
    "fortran/eftcamb/02_equispaced_interpolation_linear_1D.f90",
    "ind = int((x - self%x_initial)/self%grid_width) + 1",
    "ind = min(int((x - self%x_initial)/self%grid_width) + 1, self%num_points - 1)", "all")]  # type: ignore[list-item]


# results.f90 SetParams: test the NEW parameters (P) and only when P%EFTCAMB is allocated (this%CP%EFTCAMB may be
# unallocated on the first call: -O0 evaluates both .and. operands and dereferences it)
PATCHES["cli_p_eftcamb"] = [(
    "fortran/results.f90",
    """    if (.not. allocated(P%DarkEnergy) .and. this%CP%EFTCAMB%EFTflag == 0) then
        call GlobalError('DarkEnergy not set', error_darkenergy)
    end if""",
    """    if (.not. allocated(P%DarkEnergy)) then
        if (.not. allocated(P%EFTCAMB)) then
            call GlobalError('DarkEnergy not set', error_darkenergy)
        else if (P%EFTCAMB%EFTflag == 0) then
            call GlobalError('DarkEnergy not set', error_darkenergy)
        end if
    end if""")]


# crash-hunting continuation only: compute the diagnostic quasi-static mu/Sigma only when the cache is printed
# (they are not used by the evolution; their closed-form ratio is 0/0 when all EFT functions entering it vanish)
PATCHES["qs_printing_only"] = [(
    "fortran/equations.f90",
    """            call CP%EFTCAMB%model%compute_QS_mu(a, k, CP%eft_par_cache, EV%eft_cache)
            call CP%EFTCAMB%model%compute_QS_Sigma(a, k, CP%eft_par_cache, EV%eft_cache)
            call CP%EFTCAMB%model%compute_QS_Mu_Sigma(a, k, CP%eft_par_cache, EV%eft_cache)""",
    """            if (EV%EFTCAMB_printing) then
            call CP%EFTCAMB%model%compute_QS_mu(a, k, CP%eft_par_cache, EV%eft_cache)
            call CP%EFTCAMB%model%compute_QS_Sigma(a, k, CP%eft_par_cache, EV%eft_cache)
            call CP%EFTCAMB%model%compute_QS_Mu_Sigma(a, k, CP%eft_par_cache, EV%eft_cache)
            end if""")]


# instrumentation: report derivative components that derivs leaves unassigned (sentinel survives) or non-finite
# (also needs "use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_get_halting_mode, ieee_set_halting_mode,
#  ieee_invalid" added after "use MassiveNu" at the top of derivs)
PATCHES["nan_probe"] = [
    ("fortran/equations.f90",
     """    k = EV%k_buf
    k2 = EV%k2_buf

    ! Get background scale factor, sound speed and ionisation fraction.
    if (EV%TightCoupling) then
        call EV%ThermoData%Values(tau, a, cs2, opacity, dopacity)""",
     """    k = EV%k_buf
    k2 = EV%k2_buf
    ayprime(1:n) = -7.77e300_dl

    ! Get background scale factor, sound speed and ionisation fraction.
    if (EV%TightCoupling) then
        call EV%ThermoData%Values(tau, a, cs2, opacity, dopacity)"""),
    ("fortran/equations.f90",
     """        end if

    end if

    end subroutine derivs""",
     """        end if

    end if
    call ieee_get_halting_mode(ieee_invalid, audit_halt)
    call ieee_set_halting_mode(ieee_invalid, .false.)
    do audit_i = 1, n
        if (.not. ieee_is_finite(ayprime(audit_i)) .or. ayprime(audit_i) == -7.77e300_dl) then
            if (audit_count < 40) write(*, '(a,i5,a,i5,a,i5,a,l2,a,l2,a,l2,a,3ES14.5)') 'AUDIT_DERIVS unset/nonfinite i=', &
                audit_i, ' n=', n, ' w_ix=', EV%w_ix, ' active=', EV%EFTCAMBactive, ' tight=', EV%TightCoupling, &
                ' nonu=', EV%no_nu_multpoles, ' tau,a,k=', tau, a, k
            audit_count = audit_count + 1
        end if
    end do
    call ieee_set_halting_mode(ieee_invalid, audit_halt)

    end subroutine derivs"""),
    ("fortran/equations.f90",
     """    real(dl) :: pigddot, E2ddot
    ! EFTCAMB MOD END.

    k = EV%k_buf""",
     """    real(dl) :: pigddot, E2ddot
    ! EFTCAMB MOD END.
    integer :: audit_i
    logical :: audit_halt
    integer, save :: audit_count = 0

    k = EV%k_buf"""),
]


# remove the dead generic z/dz evaluation in the EFT RSA branch of derivs: both sub-branches below recompute z and dz,
# but with EFTCAMB_evolve_delta_phi = T the EFTeom* coefficients are never computed (EFTeomG = EFTeomU = EFTeomQ = 0)
# and the dead statement divides by zero (FP trap at equations.f90:3008-3010 for 5_hdsk_TG_1/3)
PATCHES["rsa_dead"] = [(
    "fortran/equations.f90",
    """        ! RSA

        z = 1._dl/EV%eft_cache%EFTeomG* &
            (EV%eft_cache%EFTeomQ*etak/adotoa + 0.5_dl*dgrho/(k*adotoa*(1._dl + EV%eft_cache%EFTOmegaV)) + &
            EV%eft_cache%EFTeomL/(k*CP%eft_par_cache%h0_Mpc))
        dz = 1._dl/EV%eft_cache%EFTeomU* &
            (-2._dl*adotoa*z*(1._dl + EV%eft_cache%EFTeomY - 0.5_dl*EV%eft_cache%EFTeomG/EV%eft_cache%EFTeomQ) - &
            0.5_dl*dgrho/k/(1._dl + EV%eft_cache%EFTOmegaV)/EV%eft_cache%EFTeomQ - &
            adotoa*EV%eft_cache%EFTeomL/(k*CP%eft_par_cache%h0_Mpc)/EV%eft_cache%EFTeomQ - &
            1.5_dl/k/(1._dl + EV%eft_cache%EFTOmegaV)*EV%eft_cache%EFTeomM/CP%eft_par_cache%h0_Mpc)
""",
    """        ! RSA
""")]

# nan_probe variant that keeps ayprime(EV%w_ix + 1) (derivs reads it before writing it, see output-pidotdot-zero.py)
PATCHES["nan_probe2"] = [list(t) for t in PATCHES["nan_probe"]]  # type: ignore[misc]
PATCHES["nan_probe2"][0][2] = PATCHES["nan_probe2"][0][2].replace(
    "    ayprime(1:n) = -7.77e300_dl\n",
    "    audit_keep = 0._dl\n    if (EV%w_ix > 0 .and. EV%w_ix + 1 <= n) audit_keep = ayprime(EV%w_ix + 1)\n"
    "    ayprime(1:n) = -7.77e300_dl\n    if (EV%w_ix > 0 .and. EV%w_ix + 1 <= n) ayprime(EV%w_ix + 1) = audit_keep\n")
PATCHES["nan_probe2"][2][2] = PATCHES["nan_probe2"][2][2].replace(
    "    integer :: audit_i\n", "    integer :: audit_i\n    real(dl) :: audit_keep\n")
PATCHES["nan_probe2"].append((
    "fortran/equations.f90",
    """    use constants, only: barssc0, Compton_CT, line21_const
    use MassiveNu
""",
    """    use constants, only: barssc0, Compton_CT, line21_const
    use MassiveNu
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_get_halting_mode, ieee_set_halting_mode, ieee_invalid
"""))
PATCHES["nan_probe2"] = [tuple(t) for t in PATCHES["nan_probe2"]]  # type: ignore[misc]


def apply(tree: Path, name: str) -> None:
    for item in PATCHES[name]:
        rel, old, new = item[:3]
        mode = item[3] if len(item) > 3 else "once"
        p = tree / rel
        s = p.read_text()
        n = s.count(old)
        if mode == "once" and n != 1:
            raise SystemExit(f"{name}: expected 1 match in {rel}, found {n}")
        if mode == "all" and n < 1:
            raise SystemExit(f"{name}: no match in {rel} for {old[:60]!r}")
        p.write_text(s.replace(old, new))
        print(f"{name}: {rel} ({n} replacement{'s' if n != 1 else ''})")


if __name__ == "__main__":
    tree = Path(sys.argv[1])
    for nm in sys.argv[2:]:
        apply(tree, nm)

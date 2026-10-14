! ----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2023 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
! ----------------------------------------------------------------------------------------

!> @file 11_EFTCAMB_python.f90
!! This file contains the code to interface the EFTCAMB code to the python side
!! of CAMB. There is no physics here, just wrappers...


! ----------------------------------------------------------------------------------------
!> This module contains the code to interface the EFTCAMB code to the python side
!! of CAMB. There is no physics here, just wrappers...

!> @author Marco Raveri

    module EFTCAMB_python

    use IniObjects
    use results
    use EFTCAMB_stability
    use MassiveNu

    implicit none

    type(TIniFile), private :: Ini

    contains

    ! ---------------------------------------------------------------------------------------------
    !> Helper function to strip all whitespaces
    function helper_stripspaces(string)
    character(len=*) :: string
    character(len=:), allocatable :: helper_stripspaces
    integer :: stringLen
    integer :: last, actual
    stringLen = len (string)
    last = 1
    actual = 1
    do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        end if
    end do
    helper_stripspaces = trim(string)
    end function

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to initialize EFTCAMB from python
    subroutine initialize_EFTCAMB_from_py_read_params(key, nkey, value, nvalue, mode) bind(c, &
        name='initialize_EFTCAMB_from_py_read_params_' )

    use iso_c_binding

    implicit none

    ! interface:
    integer(c_int), value, intent(in) :: nkey, nvalue, mode
    character(len=1), intent(in) :: key(nkey), value(nvalue)
    ! local fortran copies:
    character(len=nkey) :: f_key
    character(len=nvalue) :: f_value
    integer :: i

    ! initialize local copies:
    do i = 1, nkey
        f_key(i:i) = key(i)(1:1)
    end do
    do i = 1, nvalue
        f_value(i:i) = value(i)(1:1)
    end do

    ! if mode == 0 then create the INI (virtual) file
    if (mode == 0) then
        call Ini%Close()
        call Ini%Init(ignoreDuplicates=.true.)
        ! if mode == 1 append settings to the INI file
    else if (mode == 1) then
        call Ini%AddString(trim(helper_stripspaces(f_key)), trim(helper_stripspaces(f_value)))
    end if

    end subroutine initialize_EFTCAMB_from_py_read_params

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to initialize EFTCAMB from python
    subroutine initialize_EFTCAMB_from_py(State, eft_error)

    implicit none

    type(CAMBdata) :: State
    integer        :: eft_error
    logical        :: reset_in_place
    type(TEFTCAMB) :: default_EFTCAMB
    type(TEFTCAMB_parameter_cache) :: default_eft_par_cache

    ! allocate, or reset to the default values in place. Keeping the objects at the same address
    ! keeps valid the python handles to them (e.g. eftcamb = pars.EFTCAMB) across initializations.
    reset_in_place = .false.
    if (allocated(State%CP%EFTCAMB)) then
        select type (eft => State%CP%EFTCAMB)
        type is (TEFTCAMB)
            eft = default_EFTCAMB
            reset_in_place = .true.
        end select
    end if
    if (.not. reset_in_place) then
        if (allocated(State%CP%EFTCAMB)) deallocate(State%CP%EFTCAMB)
        allocate(State%CP%EFTCAMB)
    end if
    reset_in_place = .false.
    if (allocated(State%CP%eft_par_cache)) then
        select type (cache => State%CP%eft_par_cache)
        type is (TEFTCAMB_parameter_cache)
            cache = default_eft_par_cache
            reset_in_place = .true.
        end select
    end if
    if (.not. reset_in_place) then
        if (allocated(State%CP%eft_par_cache)) deallocate(State%CP%eft_par_cache)
        allocate(State%CP%eft_par_cache)
    end if
    ! read the options in:
    call State%CP%EFTCAMB%EFTCAMB_init_from_file(Ini)
    ! call the eftcamb allocation:
    if (State%CP%EFTCAMB%EFTflag /= 0) then
        ! initialize the model from file:
        call State%CP%EFTCAMB%EFTCAMB_init_model_from_file(Ini, eft_error)
    end if
    ! store the parameters that have been read in the EFTCAMB object:
    call State%CP%EFTCAMB%EFTCAMB_store_read_parameters(Ini)

    end subroutine initialize_EFTCAMB_from_py

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get the keys and values of the parameters read when an EFTCAMB object was initialized
    subroutine get_read_parameters(EFTCAMB, key, value, index, mode)

    implicit none

    type(TEFTCAMB) :: EFTCAMB
    character(len=200) :: key, value
    !< if mode=0 this returns the number of parameters, otherwise this is used as the parameter number to return
    integer :: index
    !< if mode=0 returns the number of parameters only, if mode=1 then returns the key, value combination
    integer :: mode

    ! if mode == 0 return the number of elements
    if (mode == 0) then
        index = 0
        if (allocated(EFTCAMB%read_param_keys)) index = size(EFTCAMB%read_param_keys)
        ! if mode == 1 return the key, value combination
    else if (mode == 1) then
        key = EFTCAMB%read_param_keys(index)
        value = EFTCAMB%read_param_values(index)
    end if

    end subroutine get_read_parameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to print EFTCAMB feedback:
    subroutine EFTCAMB_feedback(EFTCAMB, print_params)

    implicit none

    type(TEFTCAMB) :: EFTCAMB
    logical        :: print_params

    ! print the EFTCAMB header:
    call EFTCAMB%EFTCAMB_print_header()
    ! print model feedback:
    call EFTCAMB%EFTCAMB_print_model_feedback(print_params)

    end subroutine EFTCAMB_feedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get number of model parameters:
    subroutine EFTCAMB_get_num_params(EFTCAMB, num_params)

    implicit none

    type(TEFTCAMB) :: EFTCAMB
    integer        :: num_params

    num_params = EFTCAMB%model%parameter_number

    end subroutine EFTCAMB_get_num_params

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get parameter names:
    subroutine EFTCAMB_get_param_names(EFTCAMB, i, name)

    implicit none

    type(TEFTCAMB)     :: EFTCAMB
    integer            :: i
    character(len=200) :: name

    call EFTCAMB%model%parameter_names(i, name)

    end subroutine EFTCAMB_get_param_names

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get parameter labels:
    subroutine EFTCAMB_get_param_labels(EFTCAMB, i, name)

    implicit none

    type(TEFTCAMB)     :: EFTCAMB
    integer            :: i
    character(len=200) :: name

    call EFTCAMB%model%parameter_names_latex(i, name)

    end subroutine EFTCAMB_get_param_labels

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get parameter values:
    subroutine EFTCAMB_get_param_values(EFTCAMB, i, value)

    implicit none

    type(TEFTCAMB) :: EFTCAMB
    integer        :: i
    real(dl)       :: value

    call EFTCAMB%model%parameter_values(i, value)

    end subroutine EFTCAMB_get_param_values

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get the EFTCAMB model name:
    subroutine EFTCAMB_get_model_name(EFTCAMB, model_name)

    implicit none

    type(TEFTCAMB)     :: EFTCAMB
    character(len=200) :: model_name

    if (allocated(EFTCAMB%model)) then
        model_name = EFTCAMB%model%name
    else
        model_name = ''
    end if

    end subroutine EFTCAMB_get_model_name

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to set the EFTCAMB model name:
    subroutine EFTCAMB_set_model_name(EFTCAMB, model_name)

    implicit none

    type(TEFTCAMB)     :: EFTCAMB
    character(len=200) :: model_name

    if (allocated(EFTCAMB%model)) then
        EFTCAMB%model%name = trim(adjustl(model_name))
    end if

    end subroutine EFTCAMB_set_model_name

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get EFTCAMB effective w0wa values:
    subroutine EFTCAMB_get_effective_w0wa_values(EFTCAMB, w0, wa)

    implicit none

    type(TEFTCAMB) :: EFTCAMB
    real(dl)       :: w0        !< the effective w0 value (output)
    real(dl)       :: wa        !< the effective wa value (output)

    call EFTCAMB%Effective_w_wa(w0, wa)

    end subroutine EFTCAMB_get_effective_w0wa_values

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get the EFT time evolution:
    subroutine EFTCAMB_GetEvolution(this, nq, q, neta, eta, timestep_cache)

    use GaugeInterface
    use CAMBmain

    type(CAMBdata), target :: this
    integer, intent(in)  :: nq, neta
    real(dl), intent(in) :: q(nq), eta(neta)
    type(TEFTCAMB_timestep_cache) :: timestep_cache(neta*nq)

    real(dl) :: taustart
    integer  :: q_ix
    type(EvolutionVars) :: Ev
    ! activate status:
    call SetActiveState(this)
    ! initialize:
    global_error_flag = 0
    ! get initial time and call Thermo initialization if needed:
    taustart = min(minval(eta), GetTauStart(maxval(q)))
    if (.not. this%ThermoData%HasThermoData .or. taustart < this%ThermoData%tauminn) call this%ThermoData%Init(this, &
        taustart)
    ! do the parallel loop over the requested k modes.
    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(DYNAMIC), PRIVATE(EV, q_ix)
    do q_ix = 1, nq
        if (global_error_flag == 0) then
            ! initialize evolution variables:
            Ev%q_ix = q_ix
            Ev%q = q(q_ix)
            Ev%TransferOnly = .false.
            Ev%q2 = Ev%q**2
            Ev%ThermoData => this%ThermoData
            call GetNumEqns(Ev)
            ! call utility function:
            call GetEFTOutputEvolutionFork(State, Ev, neta, eta, timestep_cache((q_ix - 1)*neta + 1:q_ix*neta))
        end if
    end do
    !$OMP END PARALLEL DO

    end subroutine EFTCAMB_GetEvolution

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get the EFT time evolution k mode by kmode:
    subroutine GetEFTOutputEvolutionFork(this, EV, neta, eta, timestep_cache)

    use GaugeInterface
    use CAMBmain
    use RungeKuttaDP45Module, only: RungeKuttaDP45Settings

    type(CAMBdata) :: this
    type(EvolutionVars) :: EV
    integer, intent(in) :: neta
    real(dl), intent(in) :: eta(neta)
    type(TEFTCAMB_timestep_cache) :: timestep_cache(neta)

    type(RungeKuttaDP45Settings) :: rk_settings
    real(dl) w(EV%nvar, 9), y(EV%nvar), yprime(EV%nvar)
    real(dl) tau, tol1, tauend, taustart
    integer j, ind

    ! intialize:
    w = 0
    y = 0
    taustart = min(GetTauStart(min(500._dl, EV%q)), minval(eta))
    ! call initial conditions:
    call initial(EV, y, taustart)
    ! do the time loop:
    tau = taustart
    ind = 1
    tol1 = base_tol/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost - 1)
    do j = 1, neta
        tauend = eta(j)
        if (tauend < taustart) cycle
        ! evolve the preturbations:
        EV%EFTCAMB_printing = .false.
        call GaugeInterface_EvolveScal(EV, tau, y, tauend, tol1, ind, rk_settings, w)
        yprime = 0
        EV%EFTCAMB_printing = .true.
        call EV%eft_cache%initialize()
        call derivs(EV, EV%ScalEqsToPropagate, tau, y, yprime)
        ! now EV%eft_cache is filled with all quantities at the required time:
        timestep_cache(j) = EV%eft_cache
        if (global_error_flag /= 0) return
    end do

    end subroutine GetEFTOutputEvolutionFork

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get the EFT functions at scale factor a (no recomputation):
    subroutine EFTCAMB_GetEFTFunctions(State, na, as, etas, timestep_cache)
    ! Background quantities and EFT functions at the requested scale factors. The calculation follows
    ! the background part of derivs in equations.f90, so that the results are the same as get_evolution.
    type(CAMBdata)                :: State
    integer, intent(in)           :: na
    real(dl), intent(in)          :: as(na)
    real(dl), intent(in)          :: etas(na)
    type(TEFTCAMB_timestep_cache) :: timestep_cache(na)

    integer i, nu_i
    real(dl) :: grhonu, gpnu, grhonu_i, gpnu_i, grhormass_t, gpnudot_i, am
    real(dl) :: a, adotoa, grhorad

    associate(EFTCAMB => State%CP%EFTCAMB, par_cache => State%CP%eft_par_cache)
        do i = 1, na
            associate(cache => timestep_cache(i))
                a = as(i)
                cache%a = a
                cache%tau = etas(i)
                ! background densities and pressure of the matter species:
                cache%grhob_t = par_cache%grhob/a
                cache%grhoc_t = par_cache%grhoc/a
                cache%grhor_t = par_cache%grhornomass/a/a
                cache%grhog_t = par_cache%grhog/a/a
                grhorad = cache%grhog_t + cache%grhor_t
                grhonu = 0._dl
                gpnu = 0._dl
                if (par_cache%Num_Nu_Massive /= 0) then
                    do nu_i = 1, par_cache%Nu_mass_eigenstates
                        grhormass_t = par_cache%grhormass(nu_i)/a/a
                        call ThermalNuBackground%rho_P(a*par_cache%nu_masses(nu_i), grhonu_i, gpnu_i)
                        grhonu = grhonu + grhormass_t*grhonu_i
                        gpnu = gpnu + grhormass_t*gpnu_i
                    end do
                end if
                cache%grhom_t = cache%grhob_t + cache%grhoc_t + grhorad + grhonu
                cache%grhoa2 = cache%grhom_t/a**4
                cache%gpresm_t = gpnu + grhorad/3._dl
                ! expansion history:
                if (.not. EFTCAMB%EFTCAMB_model_is_designer) then
                    call EFTCAMB%model%compute_background_EFT_functions(a, par_cache, cache)
                    call EFTCAMB%model%compute_adotoa(a, par_cache, cache)
                else
                    call EFTCAMB%model%compute_adotoa(a, par_cache, cache)
                end if
                adotoa = cache%adotoa
                ! massive neutrinos and first derivative of the matter pressure:
                cache%grhonu_tot = 0._dl
                cache%gpinu_tot = 0._dl
                cache%grhonudot_tot = 0._dl
                cache%gpinudot_tot = 0._dl
                cache%gpinudotdot_tot = 0._dl
                if (par_cache%Num_Nu_Massive /= 0) then
                    do nu_i = 1, par_cache%Nu_mass_eigenstates
                        grhormass_t = par_cache%grhormass(nu_i)/a/a
                        am = a*par_cache%nu_masses(nu_i)
                        call ThermalNuBackground%rho_P(am, grhonu_i, gpnu_i)
                        cache%grhonu_tot = cache%grhonu_tot + grhormass_t*grhonu_i
                        cache%gpinu_tot = cache%gpinu_tot + grhormass_t*gpnu_i
                        cache%gpinudot_tot = cache%gpinudot_tot + &
                            grhormass_t*(ThermalNuBackground%pidot(am, adotoa, gpnu_i) - 4._dl*adotoa*gpnu_i)
                        cache%grhonudot_tot = cache%grhonudot_tot + &
                            grhormass_t*(ThermalNuBackground%drho(am, adotoa) - 4._dl*adotoa*grhonu_i)
                    end do
                end if
                cache%gpresdotm_t = -4._dl*adotoa*grhorad/3._dl + cache%gpinudot_tot
                call EFTCAMB%model%compute_H_derivs(a, par_cache, cache)
                ! second derivative of the massive neutrino and matter pressure (needs Hdot):
                if (par_cache%Num_Nu_Massive /= 0) then
                    do nu_i = 1, par_cache%Nu_mass_eigenstates
                        grhormass_t = par_cache%grhormass(nu_i)/a/a
                        am = a*par_cache%nu_masses(nu_i)
                        call ThermalNuBackground%rho_P(am, grhonu_i, gpnu_i)
                        gpnudot_i = ThermalNuBackground%pidot(am, adotoa, gpnu_i)
                        cache%gpinudotdot_tot = cache%gpinudotdot_tot - &
                            4._dl*adotoa*grhormass_t*(gpnudot_i - 4._dl*adotoa*gpnu_i) + &
                            grhormass_t*(ThermalNuBackground%pidotdot(am, adotoa, cache%Hdot, gpnu_i, gpnudot_i) &
                            - 4._dl*cache%Hdot*gpnu_i - 4._dl*adotoa*gpnudot_i)
                    end do
                end if
                cache%gpresdotdotm_t = -(4._dl/3._dl)*cache%Hdot*grhorad + (16._dl/3._dl)*grhorad*adotoa**2 + &
                    cache%gpinudotdot_tot
                ! third derivative of H:
                call EFTCAMB%model%compute_H_derivs(a, par_cache, cache)
                ! remaining background quantities:
                if (EFTCAMB%EFTCAMB_model_is_designer) then
                    call EFTCAMB%model%compute_background_EFT_functions(a, par_cache, cache)
                end if
                call EFTCAMB%model%compute_rhoQPQ(a, par_cache, cache)
                call EFTCAMB%model%compute_secondorder_EFT_functions(a, par_cache, cache)
            end associate
        end do
    end associate

    end subroutine EFTCAMB_GetEFTFunctions

    ! ---------------------------------------------------------------------------------------------

    end module EFTCAMB_python

! ----------------------------------------------------------------------------------------

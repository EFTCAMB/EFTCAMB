!----------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------

!> @file 11_EFTCAMB_python.f90
!! This file contains the code to interface the EFTCAMB code to the python side
!! of CAMB. There is no physics here, just wrappers...


!----------------------------------------------------------------------------------------
!> This module contains the code to interface the EFTCAMB code to the python side
!! of CAMB. There is no physics here, just wrappers...

!> @author Marco Raveri

module EFTCAMB_python

    use IniObjects
    use results
    use EFTCAMB_stability

    implicit none

    type(TIniFile), private   :: Ini

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
          endif
      end do
      helper_stripspaces = trim(string)
  end function

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to initialize EFTCAMB from python
  subroutine initialize_EFTCAMB_from_py_read_params( key, nkey, value, nvalue, mode ) bind( c, name='initialize_EFTCAMB_from_py_read_params_' )

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
    do i=1, nkey
      f_key(i:i) = key(i)(1:1)
    end do
    do i=1, nvalue
      f_value(i:i) = value(i)(1:1)
    end do

    ! if mode == 0 then create the INI (virtual) file
    if (mode == 0) then
        call Ini%Close()
        call Ini%Init( ignoreDuplicates=.True. )
    ! if mode == 1 append settings to the INI file
    else if (mode == 1) then
        call Ini%AddString( trim(helper_stripspaces(f_key)), trim(helper_stripspaces(f_value)) )
    end if

  end subroutine initialize_EFTCAMB_from_py_read_params

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to initialize EFTCAMB from python
  subroutine initialize_EFTCAMB_from_py( State, eft_error )

    implicit none

    Type(CAMBdata)  :: State
    integer         :: eft_error

    ! allocate:
    if ( allocated(State%CP%EFTCAMB) ) deallocate(State%CP%EFTCAMB)
    allocate( State%CP%EFTCAMB )
    if ( allocated(State%CP%eft_par_cache) ) deallocate(State%CP%eft_par_cache)
    allocate( State%CP%eft_par_cache )
    ! read the options in:
    call State%CP%EFTCAMB%EFTCAMB_init_from_file(Ini)
    ! call the eftcamb allocation:
    if ( State%CP%EFTCAMB%EFTFlag /= 0 ) then
      ! initialize the model from file:
      call State%CP%EFTCAMB%EFTCAMB_init_model_from_file(Ini, eft_error)
    end if

  end subroutine initialize_EFTCAMB_from_py

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to get the keys and values that are in the ini file
  subroutine get_read_parameters( key, value, index, mode )

    implicit none

    character(len=200) :: key, value
    integer :: index !< if mode=0 this returns the number of parameters, otherwise this is used as the parameter number to return
    integer :: mode  !< if mode=0 returns the number of parameters only, if mode=1 then returns the key, value combination

    ! if mode == 0 return the number of elements
    if (mode == 0) then
      index = Ini%ReadValues%Count
    ! if mode == 1 append settings to the INI file
    else if (mode == 1) then
      key   = Ini%ReadValues%Items(index)%P%Name
      value = Ini%ReadValues%Items(index)%P%Value
    end if

  end subroutine get_read_parameters

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to print EFTCAMB feedback:
  subroutine EFTCAMB_feedback( EFTCAMB, print_params )

    implicit none

    Type(TEFTCAMB) :: EFTCAMB
    logical        :: print_params

    ! print the EFTCAMB header:
    call EFTCAMB%EFTCAMB_print_header()
    ! print model feedback:
    call EFTCAMB%EFTCAMB_print_model_feedback(print_params)

  end subroutine EFTCAMB_feedback

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to get number of model parameters:
  subroutine EFTCAMB_get_num_params( EFTCAMB, num_params )

    implicit none

    Type(TEFTCAMB) :: EFTCAMB
    integer        :: num_params

    num_params = EFTCAMB%model%parameter_number

  end subroutine EFTCAMB_get_num_params

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to get parameter names:
  subroutine EFTCAMB_get_param_names( EFTCAMB, i, name )

    implicit none

    Type(TEFTCAMB)     :: EFTCAMB
    integer            :: i
    character(len=200) :: name

    call EFTCAMB%model%parameter_names( i, name )

  end subroutine EFTCAMB_get_param_names

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to get parameter labels:
  subroutine EFTCAMB_get_param_labels( EFTCAMB, i, name )

    implicit none

    Type(TEFTCAMB)     :: EFTCAMB
    integer            :: i
    character(len=200) :: name

    call EFTCAMB%model%parameter_names_latex( i, name )

  end subroutine EFTCAMB_get_param_labels

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to get parameter values:
  subroutine EFTCAMB_get_param_values( EFTCAMB, i, value )

    implicit none

    Type(TEFTCAMB) :: EFTCAMB
    integer        :: i
    real(dl)       :: value

    call EFTCAMB%model%parameter_values(i, value)

  end subroutine EFTCAMB_get_param_values

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to get the EFTCAMB model name:
  subroutine EFTCAMB_get_model_name( EFTCAMB, model_name )

    implicit none

    Type(TEFTCAMB)     :: EFTCAMB
    character(len=200) :: model_name

    if (allocated(EFTCAMB%model)) then
      model_name = EFTCAMB%model%name
    else
      model_name = ''
    end if

  end subroutine EFTCAMB_get_model_name

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to set the EFTCAMB model name:
  subroutine EFTCAMB_set_model_name( EFTCAMB, model_name )

    implicit none

    Type(TEFTCAMB)     :: EFTCAMB
    character(len=200) :: model_name

    if (allocated(EFTCAMB%model)) then
      EFTCAMB%model%name = trim(adjustl(model_name))
    end if

  end subroutine EFTCAMB_set_model_name

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to get EFTCAMB effective w0wa values:
  subroutine EFTCAMB_get_effective_w0wa_values( EFTCAMB, w0, wa )

    implicit none

    Type(TEFTCAMB)     :: EFTCAMB
    real(dl)           :: w0        !< the effective w0 value (output)
    real(dl)           :: wa        !< the effective wa value (output)

    call EFTCAMB%Effective_w_wa(w0, wa)
    
  end subroutine EFTCAMB_get_effective_w0wa_values

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to get the EFT time evolution:
  subroutine EFTCAMB_GetEvolution( this, nq, q, neta, eta, timestep_cache )

    use GaugeInterface
    use CAMBmain

    Type(CAMBdata),target :: this
    integer, intent(in)  :: nq, neta
    real(dl), intent(in) :: q(nq), eta(neta)
    type(TEFTCAMB_timestep_cache) :: timestep_cache(neta*nq)

    real(dl) :: taustart
    integer  :: q_ix
    Type(EvolutionVars) :: Ev
    ! activate status:
    call SetActiveState(this)
    ! initialize:
    global_error_flag = 0
    ! get initial time and call Thermo initialization if needed:
    taustart = min(minval(eta),GetTauStart(maxval(q)))
    if (.not. this%ThermoData%HasThermoData .or. taustart < this%ThermoData%tauminn) call this%ThermoData%Init(this,taustart)
    ! do the parallel loop over the requested k modes.
    !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC), PRIVATE(EV, q_ix)
    do q_ix = 1, nq
        if (global_error_flag==0) then
            ! initialize evolution variables:
            EV%q_ix = q_ix
            EV%q = q(q_ix)
            EV%TransferOnly=.false.
            EV%q2=EV%q**2
            EV%ThermoData => this%ThermoData
            call GetNumEqns(EV)
            ! call utility function:
            call GetEFTOutputEvolutionFork( State, EV, neta, eta, timestep_cache((q_ix-1)*neta+1:q_ix*neta) )
        end if
    end do
    !$OMP END PARALLEL DO

  end subroutine EFTCAMB_GetEvolution

  ! ---------------------------------------------------------------------------------------------
  !> Subroutine to get the EFT time evolution k mode by kmode:
  subroutine GetEFTOutputEvolutionFork(this, EV, neta, eta, timestep_cache)

    use GaugeInterface
    use CAMBmain

    type(CAMBdata) :: this
    type(EvolutionVars) :: EV
    integer, intent(in)  :: neta
    real(dl), intent(in) :: eta(neta)
    type(TEFTCAMB_timestep_cache) :: timestep_cache(neta)

    real(dl) c(24), w(EV%nvar,9), y(EV%nvar), yprime(EV%nvar)
    real(dl) tau, tol1, tauend, taustart
    integer j,ind

    ! intialize:
    w=0
    y=0
    taustart = min(GetTauStart(min(500._dl,EV%q)),minval(eta))
    ! call initial conditions:
    call initial(EV,y, taustart)
    ! do the time loop:
    tau=taustart
    ind=1
    tol1=tol/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost-1)
    do j = 1, neta
        tauend = eta(j)
        if (tauend<taustart) cycle
        ! evolve the preturbations:
        EV%EFTCAMB_printing = .False.
        call GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
        yprime = 0
        EV%EFTCAMB_printing = .True.
        call EV%eft_cache%initialize()
        call derivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)
        ! now EV%eft_cache is filled with all quantities at the required time:
        timestep_cache(j) = EV%eft_cache
        if (global_error_flag/=0) return
    end do

  end subroutine GetEFTOutputEvolutionFork

  ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_python

!----------------------------------------------------------------------------------------

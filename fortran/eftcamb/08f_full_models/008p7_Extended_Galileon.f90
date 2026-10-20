! ----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2017 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
! ----------------------------------------------------------------------------------------

!> @file 10p5_Extended_Galileon.f90
!! This file contains the definition of the Extended Galileon model.
!! Please refer to the numerical notes for details.

! ----------------------------------------------------------------------------------------
!> This module contains the definition of the Extended Galileon model.
!! Please refer to the numerical notes for details.

!> @author Luis Atayde, Noemi Frusciante, Simone Peirone

!> updated to pyEFTCAMB by Fabrizio Renzi

    module EFTCAMB_full_Extended_Galileon

    use precision
    use config, only: GlobalError, error_unsupported_params
    use IniObjects
    use MpiUtils
    use FileUtils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_model_full
    use MassiveNu

    implicit none

    private

    public EFTCAMB_Extended_Galileon

    ! ----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of the Extended Galileon model.
    type, extends (EFTCAMB_full_model) :: EFTCAMB_Extended_Galileon

        ! the model parameters:
        real(dl) :: ExtendedGalileon_B      !< Extended Galileon model parameter \f$B\f$
        real(dl) :: ExtendedGalileon_q      !< Extended Galileon model parameter \f$q\f$
        !< Extended Galileon background parameter \f$\xi\f$ deriving from the tracker solution
        real(dl) :: csi
        real(dl) :: S

    contains

    ! initialization of the model:
    !< subroutine that reads the parameters of the model from file
    procedure :: read_model_selection => EFTCAMBExtendedGalileonReadModelSelectionFromFile
    !< subroutine that allocates the model selection. For Horava this is a dummy procedure.
    procedure :: allocate_model_selection => EFTCAMBExtendedGalileonAllocateModelSelection
    !< subroutine that initializes the model parameters based on the values found in an input array.
    procedure :: init_model_parameters => EFTCAMBExtendedGalileonInitModelParameters
    !< subroutine that reads the parameters of the model from file.
    procedure :: init_model_parameters_from_file => EFTCAMBExtendedGalileonInitModelParametersFromFile

    ! background solver:
    !< subroutine that initializes the background of Extended Galileon.
    procedure :: initialize_background => EFTCAMBExtendedGalileonInitBackground

    ! utility functions:
    !< subroutine that computes the number of parameters of the model.
    procedure :: compute_param_number => EFTCAMBExtendedGalileonComputeParametersNumber
    !< subroutine that prints on the screen feedback information about the model.
    procedure :: feedback => EFTCAMBExtendedGalileonFeedback
    !< subroutine that returns the i-th parameter name of the model.
    procedure :: parameter_names => EFTCAMBExtendedGalileonParameterNames
    !< subroutine that returns the i-th parameter name of the model.
    procedure :: parameter_names_latex => EFTCAMBExtendedGalileonParameterNamesLatex
    !< subroutine that returns the i-th parameter value.
    procedure :: parameter_values => EFTCAMBExtendedGalileonParameterValues

    ! CAMB related procedures:
    !< subroutine that computes the value of the background EFT functions at a given time.
    procedure :: compute_background_EFT_functions => EFTCAMBExtendedGalileonBackgroundEFTFunctions
    !< subroutine that computes the value of the second order EFT functions at a given time.
    procedure :: compute_secondorder_EFT_functions => EFTCAMBExtendedGalileonSecondOrderEFTFunctions
    !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    procedure :: compute_adotoa => EFTCAMBExtendedGalileonComputeAdotoa
    !< subroutine that computes the two derivatives wrt conformal time of H.
    procedure :: compute_H_derivs => EFTCAMBExtendedGalileonComputeHubbleDer

    ! stability procedures:
    !< function that computes model specific stability requirements.
    procedure :: additional_model_stability => EFTCAMBExtendedGalileonAdditionalModelStability

    end type EFTCAMB_Extended_Galileon

    ! ---------------------------------------------------------------------------------------------

    contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBExtendedGalileonReadModelSelectionFromFile(self, Ini, eft_error)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self          !< the base class
    class(TIniFile)                  :: Ini           !< Input ini file
    integer                          :: eft_error     !< error code: 0 all fine, 1 initialization failed

    end subroutine EFTCAMBExtendedGalileonReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBExtendedGalileonAllocateModelSelection(self, Ini, eft_error)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self        !< the base class
    class(TIniFile)                  :: Ini         !< Input ini file
    integer                          :: eft_error   !< error code: 0 all fine, 1 initialization failed

    ! Nothing needs to be done but procedure present because its definition is deferred.

    end subroutine EFTCAMBExtendedGalileonAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBExtendedGalileonInitModelParameters(self, array)

    implicit none

    class(EFTCAMB_Extended_Galileon)                       :: self   !< the base class
    !< input array with the values of the parameters of the model.
    real(dl), dimension(self%parameter_number), intent(in) :: array

    self%S = array(1)
    self%ExtendedGalileon_q = array(2)

    end subroutine EFTCAMBExtendedGalileonInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBExtendedGalileonInitModelParametersFromFile(self, Ini, eft_error)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self      !< the base class
    class(TIniFile)                  :: Ini       !< Input ini file
    integer                          :: eft_error !< error code: 0 all fine, 1 initialization failed

    self%S = Ini%Read_Double('ExtendedGalileon_S', 0._dl)
    self%ExtendedGalileon_q = Ini%Read_Double('ExtendedGalileon_q', 0._dl)

    end subroutine EFTCAMBExtendedGalileonInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBExtendedGalileonComputeParametersNumber(self)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self               !< the base class

    self%parameter_number = 2        ! luis-change::0->2

    end subroutine EFTCAMBExtendedGalileonComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBExtendedGalileonFeedback(self, print_params)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self         !< the base class
    logical, optional                :: print_params !< optional flag that decised whether to print numerical values
    !! of the parameters.

    ! print general model informations:
    if (self%ExtendedGalileon_B /= 0._dl) then

        write(*, *)
        write(*, '(a,a)')       '   Model              =', self%name
        write(*, '(a,I3)')      '   Number of params   =', self%parameter_number
        write(*, '(a,F12.6)')   '                xi    =', self%csi
        write(*, '(a,F12.6)')   '                B     =', self%ExtendedGalileon_B
        write(*, '(a,F12.6)')   '                q     =', self%ExtendedGalileon_q
        write(*, '(a,F12.6)')   '                S     =', self%S

    end if

    end subroutine EFTCAMBExtendedGalileonFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBExtendedGalileonParameterNames(self, i, name)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self   !< the base class
    integer, intent(in) :: i      !< the index of the parameter
    character(*), intent(out) :: name   !< the output name of the i-th parameter

    ! check the input index:
    if (i > self%parameter_number) then
        write(*, *) 'Illegal index for parameter_names.'
        write(*, *) 'Maximum value is:', self%parameter_number
        call GlobalError('EFTCAMB error: illegal parameter index in EFTCAMBExtendedGalileonParameterNames', &
            error_unsupported_params)
        name = ''
        return
    end if
    ! return the appropriate name:
    if (i == 1) then
        name = 'ExtendedGalileon_S'
        return
    end if
    if (i == 2) then
        name = 'ExtendedGalileon_q'
        return
    end if
    if (i == 0) then
        name = 'noname'
        return
    end if

    end subroutine EFTCAMBExtendedGalileonParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBExtendedGalileonParameterNamesLatex(self, i, latexname)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self       !< the base class
    integer, intent(in)       :: i          !< The index of the parameter
    character(*), intent(out) :: latexname  !< the output latex name of the i-th parameter

    ! check the input index:
    if (i > self%parameter_number) then
        write(*, *) 'Illegal index for parameter_names_latex.'
        write(*, *) 'Maximum value is:', self%parameter_number
        call GlobalError('EFTCAMB error: illegal parameter index in EFTCAMBExtendedGalileonParameterNamesLatex', &
            error_unsupported_params)
        latexname = ''
        return
    end if
    ! return the appropriate name:
    if (i == 1) then
        latexname = 'S'
        return
    end if
    if (i == 2) then
        latexname = 'q'
        return
    end if
    if (i == 0) then
        latexname = 'noname'
        return
    end if

    end subroutine EFTCAMBExtendedGalileonParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBExtendedGalileonParameterValues(self, i, value)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self   !< the base class
    integer, intent(in)   :: i      !< The index of the parameter
    real(dl), intent(out) :: value  !< the output value of the i-th parameter

    ! check the input index:
    if (i > self%parameter_number) then
        write(*, *) 'Illegal index for parameter_value.'
        write(*, *) 'Maximum value is:', self%parameter_number
        call GlobalError('EFTCAMB error: illegal parameter index in EFTCAMBExtendedGalileonParameterValues', &
            error_unsupported_params)
        value = 0._dl
        return
    end if
    ! return the appropriate name:
    if (i == 0) then
        value = 0._dl
        return
    end if
    if (i == 1) then
        value = self%S
        return
    end if
    if (i == 2) then
        value = self%ExtendedGalileon_q
        return
    end if

    end subroutine EFTCAMBExtendedGalileonParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of Extended Galileon.
    subroutine EFTCAMBExtendedGalileonInitBackground(self, params_cache, feedback_level, success, outroot)

    implicit none

    class(EFTCAMB_Extended_Galileon)              :: self           !< the base class
    !< a EFTCAMB parameter cache containing cosmological parameters
    type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache
    integer, intent(in)  :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
    logical, intent(out) :: success        !< wether the background initialization succeded or not
    character(len=*), intent(in), optional :: outroot        !< the output root for the debug files

    real(dl) :: Omega_phi0

    Omega_phi0 = params_cache%omegav

    self%ExtendedGalileon_B = self%S*self%ExtendedGalileon_q   !< B=S*q!

    call self%feedback()

    ! the EFT functions contain B/q:
    success = (self%ExtendedGalileon_q /= 0._dl)
    if (.not. success .and. feedback_level > 0) write(*, *) 'EFTCAMB Extended Galileon: q = 0 is not allowed.'

    end subroutine EFTCAMBExtendedGalileonInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBExtendedGalileonBackgroundEFTFunctions(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_Extended_Galileon)              :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout)  :: eft_cache

    real(dl) :: a2, Omega_phi0

    Omega_phi0 = eft_par_cache%omegav
    a2 = a*a

    if (a <= 0._dl) return
    ! expansion history of the model's own background (tracker solution), independent of the cache content:
    call self%compute_adotoa(a, eft_par_cache, eft_cache)
    call self%compute_H_derivs(a, eft_par_cache, eft_cache)
    if (eft_cache%adotoa <= 0._dl) return
    ! compute the background EFT functions:
    eft_cache%EFTOmegaV = 0._dl
    eft_cache%EFTOmegaP = 0._dl
    eft_cache%EFTOmegaPP = 0._dl
    eft_cache%EFTOmegaPPP = 0._dl
    eft_cache%EFTc = -(a**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*self%ExtendedGalileon_B &
        *eft_cache%adotoa**(-2 - self%ExtendedGalileon_B/self%ExtendedGalileon_q) &
        *eft_par_cache%h0_Mpc**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q) &
        *Omega_phi0*(eft_cache%adotoa**2 - eft_cache%Hdot))/(2._dl*self%ExtendedGalileon_q)

    eft_cache%EFTLambda = -((a**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*eft_cache%adotoa** &
        (-2 - self%ExtendedGalileon_B/self%ExtendedGalileon_q)*eft_par_cache%h0_Mpc** &
        (2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*Omega_phi0* &
        ((3*self%ExtendedGalileon_q + self%ExtendedGalileon_B)*eft_cache%adotoa**2 - &
        self%ExtendedGalileon_B*eft_cache%Hdot))/self%ExtendedGalileon_q)

    eft_cache%EFTcdot = -(a**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*self%ExtendedGalileon_B* &
        eft_cache%adotoa**(-3 - self%ExtendedGalileon_B/self%ExtendedGalileon_q)*eft_par_cache%h0_Mpc** &
        (2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*Omega_phi0* &
        (self%ExtendedGalileon_B*eft_cache%adotoa**4 - 2*self%ExtendedGalileon_B*eft_cache%adotoa**2*eft_cache%Hdot + &
        (2*self%ExtendedGalileon_q + self%ExtendedGalileon_B)*eft_cache%Hdot**2 - &
        self%ExtendedGalileon_q*eft_cache%adotoa*eft_cache%Hdotdot))/(2._dl*self%ExtendedGalileon_q**2)

    eft_cache%EFTLambdadot = -((a**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*self%ExtendedGalileon_B* &
        eft_cache%adotoa**(-3 - self%ExtendedGalileon_B/self%ExtendedGalileon_q)*eft_par_cache%h0_Mpc** &
        (2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*Omega_phi0* &
        ((3*self%ExtendedGalileon_q + self%ExtendedGalileon_B)*eft_cache%adotoa**4 - &
        (3*self%ExtendedGalileon_q + 2*self%ExtendedGalileon_B)*eft_cache%adotoa**2*eft_cache%Hdot + &
        (2*self%ExtendedGalileon_q + self%ExtendedGalileon_B)*eft_cache%Hdot**2 - &
        self%ExtendedGalileon_q*eft_cache%adotoa*eft_cache%Hdotdot))/self%ExtendedGalileon_q**2)

    end subroutine EFTCAMBExtendedGalileonBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the second order EFT functions at a given time.
    subroutine EFTCAMBExtendedGalileonSecondOrderEFTFunctions(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self          !< the base class
    real(dl), intent(in)             :: a             !< the input scale factor
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: a2, Omega_phi0

    Omega_phi0 = eft_par_cache%omegav
    a2 = a*a

    if (a <= 0._dl) return
    ! expansion history of the model's own background (tracker solution), independent of the cache content:
    call self%compute_adotoa(a, eft_par_cache, eft_cache)
    call self%compute_H_derivs(a, eft_par_cache, eft_cache)
    if (eft_cache%adotoa <= 0._dl) return

    ! compute the second order EFT functions:
    eft_cache%EFTGamma1V = (a**(self%ExtendedGalileon_B/self%ExtendedGalileon_q)*self%ExtendedGalileon_B* &
        eft_cache%adotoa**(-2 - self%ExtendedGalileon_B/self%ExtendedGalileon_q)*eft_par_cache%h0_Mpc** &
        (self%ExtendedGalileon_B/self%ExtendedGalileon_q)*Omega_phi0* &
        ((1 + 12*self%ExtendedGalileon_q**2)*eft_cache%adotoa**2 - eft_cache%Hdot))/(4._dl*self%ExtendedGalileon_q)

    eft_cache%EFTGamma1P = (a**(-1 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*self%ExtendedGalileon_B* &
        eft_cache%adotoa**(-4 - self%ExtendedGalileon_B/self%ExtendedGalileon_q)*eft_par_cache%h0_Mpc** &
        (self%ExtendedGalileon_B/self%ExtendedGalileon_q)*Omega_phi0* &
        ((self%ExtendedGalileon_B + 12*self%ExtendedGalileon_q**2*self%ExtendedGalileon_B)*eft_cache%adotoa**4 - &
        2*(1 + 6*self%ExtendedGalileon_q**2)*self%ExtendedGalileon_B*eft_cache%adotoa**2*eft_cache%Hdot + &
        (2*self%ExtendedGalileon_q + self%ExtendedGalileon_B)*eft_cache%Hdot**2 - &
        self%ExtendedGalileon_q*eft_cache%adotoa*eft_cache%Hdotdot))/(4._dl*self%ExtendedGalileon_q**2)

    eft_cache%EFTGamma2V = -2*self%ExtendedGalileon_B*(a*eft_par_cache%h0_Mpc/eft_cache%adotoa)** &
        ((self%ExtendedGalileon_B + self%ExtendedGalileon_q)/self%ExtendedGalileon_q)*Omega_phi0

    eft_cache%EFTGamma2P = -2/self%ExtendedGalileon_q*self%ExtendedGalileon_B* &
        (self%ExtendedGalileon_B + self%ExtendedGalileon_q)*a**(self%ExtendedGalileon_B/self%ExtendedGalileon_q)* &
        eft_par_cache%h0_Mpc**((self%ExtendedGalileon_B + self%ExtendedGalileon_q)/self%ExtendedGalileon_q)* &
        eft_cache%adotoa**(-3 - (self%ExtendedGalileon_B)/self%ExtendedGalileon_q)*Omega_phi0* &
        (eft_cache%adotoa**2 - eft_cache%Hdot)

    eft_cache%EFTGamma3V = 0._dl
    eft_cache%EFTGamma3P = 0._dl
    eft_cache%EFTGamma4V = 0._dl
    eft_cache%EFTGamma4P = 0._dl
    eft_cache%EFTGamma4PP = 0._dl
    eft_cache%EFTGamma5V = 0._dl
    eft_cache%EFTGamma5P = 0._dl
    eft_cache%EFTGamma6V = 0._dl
    eft_cache%EFTGamma6P = 0._dl

    end subroutine EFTCAMBExtendedGalileonSecondOrderEFTFunctions


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H.
    subroutine EFTCAMBExtendedGalileonComputeAdotoa(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self          !< the base class
    real(dl), intent(in)             :: a             !< the input scale factor
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: up, upp, uppp

    call ExtGalTracker(self, a, eft_par_cache, eft_cache%adotoa, up, upp, uppp)

    end subroutine EFTCAMBExtendedGalileonComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the three derivatives wrt conformal time of H.
    subroutine EFTCAMBExtendedGalileonComputeHubbleDer(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self          !< the base class
    real(dl), intent(in)             :: a             !< the input scale factor
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: adotoa, up, upp, uppp

    ! with u = ln E: d adotoa/dN = adotoa (1 + u') and d/dtau = adotoa d/dN
    call ExtGalTracker(self, a, eft_par_cache, adotoa, up, upp, uppp)
    eft_cache%Hdot = adotoa**2*(1._dl + up)
    eft_cache%Hdotdot = adotoa**3*(2._dl*(1._dl + up)**2 + upp)
    eft_cache%Hdotdotdot = adotoa**4*(uppp + 7._dl*(1._dl + up)*upp + 6._dl*(1._dl + up)**3)

    end subroutine EFTCAMBExtendedGalileonComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that solves the tracker Friedmann equation of the model (arXiv:1912.07586 eq. 16)
    !! \f$ E^2 = \Omega_{\rm tot}(a) + \Omega_{\phi,0} E^{-S} \f$, with \f$ E = H/H_0 \f$, and returns adotoa and
    !! the first three derivatives of u = ln E with respect to N = ln a. \f$ \Omega_{\rm tot} \f$ contains matter,
    !! radiation and massive neutrinos, computed here.
    subroutine ExtGalTracker(self, a, eft_par_cache, adotoa, up, upp, uppp)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self          !< the base class
    real(dl), intent(in)             :: a             !< the input scale factor
    type(TEFTCAMB_parameter_cache), intent(in) :: eft_par_cache !< the EFTCAMB parameter cache
    real(dl), intent(out)            :: adotoa        !< conformal Hubble rate
    real(dl), intent(out)            :: up, upp, uppp !< derivatives of ln(H/H0) wrt ln a

    real(dl) :: S, Om_phi, Om, Om_p, Om_pp, Om_ppp, rho, rho_p, rho_pp, rho_ppp, p, p_p, p_pp
    real(dl) :: rhonu, pnu, dpnu, ddpnu, am, Om_m, Om_r, u, u_lo, u_hi, g, dg, du, W2, V, den
    integer  :: nu_i, iter

    adotoa = 0._dl
    up = 0._dl
    upp = 0._dl
    uppp = 0._dl
    if (a <= 0._dl) return
    S = self%S
    Om_phi = eft_par_cache%omegav

    ! 1) Omega_tot(N) = rho/(3 H0^2) and its derivatives wrt N. Massive neutrinos: rho' = -3 (rho + p).
    rho = 0._dl
    p = 0._dl
    p_p = 0._dl
    p_pp = 0._dl
    if (eft_par_cache%Num_Nu_Massive /= 0) then
        do nu_i = 1, eft_par_cache%Nu_mass_eigenstates
            am = a*eft_par_cache%nu_masses(nu_i)
            rhonu = 0._dl
            pnu = 0._dl
            call ThermalNuBackground%rho_P(am, rhonu, pnu)
            dpnu = ThermalNuBackground%pidot(am, 1._dl, pnu)
            ddpnu = ThermalNuBackground%pidotdot(am, 1._dl, 0._dl, pnu, dpnu)
            rho = rho + eft_par_cache%grhormass(nu_i)*rhonu/a**4
            p = p + eft_par_cache%grhormass(nu_i)*pnu/a**4
            p_p = p_p + eft_par_cache%grhormass(nu_i)*(dpnu - 4._dl*pnu)/a**4
            p_pp = p_pp + eft_par_cache%grhormass(nu_i)*(ddpnu - 8._dl*dpnu + 16._dl*pnu)/a**4
        end do
    end if
    rho_p = -3._dl*(rho + p)
    rho_pp = -3._dl*(rho_p + p_p)
    rho_ppp = -3._dl*(rho_pp + p_pp)
    Om_m = (eft_par_cache%omegac + eft_par_cache%omegab)/a**3
    Om_r = (eft_par_cache%omegag + eft_par_cache%omegar)/a**4
    Om = Om_m + Om_r + rho/(3._dl*eft_par_cache%h0_Mpc**2)
    Om_p = -3._dl*Om_m - 4._dl*Om_r + rho_p/(3._dl*eft_par_cache%h0_Mpc**2)
    Om_pp = 9._dl*Om_m + 16._dl*Om_r + rho_pp/(3._dl*eft_par_cache%h0_Mpc**2)
    Om_ppp = -27._dl*Om_m - 64._dl*Om_r + rho_ppp/(3._dl*eft_par_cache%h0_Mpc**2)

    ! 2) solve g(u) = exp(2u) - Om - Om_phi exp(-S u) = 0 for u = ln E (safeguarded Newton):
    u_lo = 0.5_dl*log(Om)
    u_hi = 0.5_dl*log(Om + Om_phi*exp(-S*u_lo))
    do iter = 1, 200
        if (exp(2._dl*u_hi) - Om - Om_phi*exp(-S*u_hi) >= 0._dl) exit
        u_hi = u_hi + 1._dl
    end do
    u = u_hi
    do iter = 1, 200
        g = exp(2._dl*u) - Om - Om_phi*exp(-S*u)
        if (g > 0._dl) then
            u_hi = u
        else
            u_lo = u
        end if
        dg = 2._dl*exp(2._dl*u) + S*Om_phi*exp(-S*u)
        du = -g/dg
        if (u + du <= u_lo .or. u + du >= u_hi .or. dg <= 0._dl) du = 0.5_dl*(u_lo + u_hi) - u
        u = u + du
        if (abs(du) < 1.e-15_dl*max(1._dl, abs(u)) .or. u_hi - u_lo < 1.e-15_dl*max(1._dl, abs(u))) exit
    end do
    adotoa = a*eft_par_cache%h0_Mpc*exp(u)

    ! 3) derivatives of u from W - V = Om with W = E^2, V = Om_phi E^(-S):
    W2 = exp(2._dl*u)
    V = Om_phi*exp(-S*u)
    den = 2._dl*W2 + S*V
    up = Om_p/den
    upp = (Om_pp - 4._dl*up**2*W2 + S**2*up**2*V)/den
    uppp = (Om_ppp - 12._dl*up*upp*W2 - 8._dl*up**3*W2 + 3._dl*S**2*up*upp*V - S**3*up**3*V)/den

    end subroutine ExtGalTracker

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBExtendedGalileonAdditionalModelStability(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_Extended_Galileon) :: self          !< the base class
    real(dl), intent(in)             :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.
    logical :: EFTCAMBExtendedGalileonAdditionalModelStability

    EFTCAMBExtendedGalileonAdditionalModelStability = .true.

    ! mass condition
    if ((self%S + 0.02_dl)*self%ExtendedGalileon_q < 0.5_dl) EFTCAMBExtendedGalileonAdditionalModelStability = .false.

    ! Strong Coupling
    if ((self%ExtendedGalileon_q*self%S - 1.)/(2.*self%ExtendedGalileon_q*(1. + self%S/2.)) > &
        0.) EFTCAMBExtendedGalileonAdditionalModelStability = .false.

    end function EFTCAMBExtendedGalileonAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

    end module EFTCAMB_full_Extended_Galileon

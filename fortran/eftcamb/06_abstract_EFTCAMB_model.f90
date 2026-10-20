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

!> @file 06_abstract_EFTCAMB_model.f90
!! This file contains the abstract definition of all the places where EFTCAMB interacts
!! with CAMB. All EFTCAMB models should inherit from this class or the two derived
!! classes contained in 06p2_abstract_EFTCAMB_full.f90 or 06p3_abstract_EFTCAMB_designer.f90


! ----------------------------------------------------------------------------------------
!> This module contains the abstract definition of all the places where EFTCAMB interacts
!! with CAMB. All EFTCAMB models should inherit from this class or the two derived
!! classes contained in 06p2_abstract_EFTCAMB_full.f90 or 06p3_abstract_EFTCAMB_designer.f90

!> @author Bin Hu, Marco Raveri, Simone Peirone

    module EFTCAMB_abstract_model

    use precision
    use IniObjects
    use MpiUtils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_mixed_algorithms

    implicit none

    private

    public EFTCAMB_model, EFTCAMBStabilityScales

    ! ----------------------------------------------------------------------------------------
    !> This is the abstract type for EFTCAMB models. As a rule, when there is a
    !! new model it should be declared as a class inheriting from EFTCAMB_model.
    !! This guarantees maximum performances as well as maximum flexibility.
    type, abstract :: EFTCAMB_model

        integer                       :: parameter_number !< number of parameters of the model.
        character(len=:), allocatable :: name             !< name of the model.
        character(len=:), allocatable :: name_latex       !< latex name of the model.

        !< flag that tells whether the model has an effective w0wa model limit.
        logical                       :: effective_w0wa = .false.
        real(dl)                      :: w0 = -1._dl
        real(dl)                      :: wa = 0._dl

    contains

    ! initialization of the model:
    !< subroutine that initializes the name and latex name of the model.
    procedure :: init => EFTCAMBModelInitialize
    !< subroutine that reads the parameters of the model from file.
    procedure(EFTCAMBModelReadModelSelectionFromFile), deferred :: read_model_selection
    !< subroutine that allocates the model selection.
    procedure(EFTCAMBModelAllocateModelSelection), deferred :: allocate_model_selection
    !< subroutine taht initializes the model parameters based on the values found in an input array.
    procedure(EFTCAMBModelInitModelParameters), deferred :: init_model_parameters
    !< subroutine that reads the parameters of the model from file.
    procedure(EFTCAMBModelInitModelParametersFromFile), deferred :: init_model_parameters_from_file

    ! utility functions:
    !< subroutine that computes the number of parameters of the model.
    procedure(EFTCAMBModelComputeParametersNumber), deferred :: compute_param_number
    !< subroutine that prints on the screen feedback information about the model.
    procedure(EFTCAMBModelFeedback), deferred :: feedback
    !< subroutine that returns the i-th parameter name of the model.
    procedure(EFTCAMBModelParameterNames), deferred :: parameter_names
    !< subroutine that returns the i-th parameter name of the model.
    procedure(EFTCAMBModelParameterNamesLatex), deferred :: parameter_names_latex
    !< subroutine that returns the i-th parameter value.
    procedure(EFTCAMBModelParameterValues), deferred :: parameter_values

    ! background initialization functions:
    !< subroutine that initializes the background of the model, if needed.
    procedure :: initialize_background => EFTCAMBModelInitBackground

    ! CAMB related procedures:
    !< subroutine that computes the value of the background EFT functions at a given time.
    procedure(EFTCAMBModelBackgroundEFTFunctions), deferred :: compute_background_EFT_functions
    !< subroutine that computes the value of the second order EFT functions at a given time.
    procedure(EFTCAMBModelSecondOrderEFTFunctions), deferred :: compute_secondorder_EFT_functions
    !< function that computes dtauda = 1/sqrt(a^2H^2).
    procedure(EFTCAMBModelComputeDtauda), deferred :: compute_dtauda
    !< subroutine that computes adotoa = H.
    procedure(EFTCAMBModelComputeAdotoa), deferred :: compute_adotoa
    !< subroutine that computes the two derivatives wrt conformal time of H.
    procedure(EFTCAMBModelComputeHubbleDer), deferred :: compute_H_derivs
    !< subroutine that computes additional derivatives of the EFT functions, as needed by the stability calculation.
    procedure(EFTCAMBModelComputeAdditionalDer), deferred :: compute_additional_derivs
    !< function that computes the conformal time difference between two scale factor
    procedure :: compute_DeltaTau => EFTCAMBModelComputeDeltaTau
    !< function that computes the cosmic time difference between two scale factor
    procedure :: compute_DeltaCosmicT => EFTCAMBModelComputeDeltaCosmicT
    !< function that computes r_s at given scale factor
    procedure :: compute_SoundHorizon => EFTCAMBModelComputeSoundHorizon

    !< subroutine that computes \rho_Q and P_Q. For details refer to the numerical notes.
    procedure :: compute_rhoQPQ => EFTCAMBModelComputeRhoQPQ
    !< subroutine that computes the Einstein equations factors. For details refer to the numerical notes.
    procedure :: compute_Einstein_factors => EFTCAMBModelComputeEinsteinFactors
    !< subroutine that computes the pi field equations factors. For details refer to the numerical notes.
    procedure :: compute_pi_factors => EFTCAMBModelComputePiFactors
    !< subroutine that computes the factors for the tensor propagation equation. For details refer to the numerical notes.
    procedure :: compute_tensor_factors => EFTCAMBModelComputeTensorFactors
    !< subroutine that computes the kinetic and gradient terms. For details refer to the numerical notes.
    procedure :: compute_stability_factors => EFTCAMBModelComputeStabilityFactors
    !< subroutine that computes the kinetic and gradient terms. For details refer to the numerical notes.
    procedure :: compute_mass_factors => EFTCAMBModelComputeMassFactors
    ! compute Sigma and mu
    !< function that computes the quasi-static expression of mu.
    procedure :: compute_QS_mu => EFTCAMBModelComputeQSMu
    !< function that computes the quasi-static expression of Sigma.
    procedure :: compute_QS_Sigma => EFTCAMBModelComputeQSSigma
    !< function that computes the quasi-static expression of Sigma.
    procedure :: compute_QS_Mu_Sigma => EFTCAMBModelComputeQSMuSigma

    ! stability procedures:
    !< function that computes model specific stability requirements.
    procedure :: additional_model_stability => EFTCAMBModelAdditionalModelStability
    ! positivity bounds:
    !< function that computes the positivity bounds. See thesis Dani de Boe for details.
    procedure :: compute_positivity_bounds => EFTCAMBModelPositivityBounds

    end type EFTCAMB_model

    ! ---------------------------------------------------------------------------------------------
    ! EFTCAMB abstract interfaces: these are all the model procedures that the user HAS to override
    ! when writing its own model.
    ! ---------------------------------------------------------------------------------------------

    abstract interface

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBModelComputeParametersNumber(self)
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model) :: self   !< the base class
    end subroutine EFTCAMBModelComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBModelReadModelSelectionFromFile(self, Ini, eft_error)
    use IniObjects
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model) :: self      !< the base class
    type(TIniFile)       :: Ini       !< Input ini file
    integer              :: eft_error !< error code: 0 all fine, 1 initialization failed
    end subroutine EFTCAMBModelReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBModelAllocateModelSelection(self, Ini, eft_error)
    use IniObjects
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model) :: self      !< the base class
    type(TIniFile)       :: Ini       !< Input ini file for additional initialization
    integer              :: eft_error !< error code: 0 all fine, 1 initialization failed
    end subroutine EFTCAMBModelAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBModelInitModelParameters(self, array)
    use precision
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model)                                   :: self   !< the base class
    real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.
    end subroutine EFTCAMBModelInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBModelInitModelParametersFromFile(self, Ini, eft_error)
    use IniObjects
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model) :: self      !< the base class
    type(TIniFile)       :: Ini       !< Input ini file
    integer              :: eft_error !< error code: 0 all fine, 1 initialization failed
    end subroutine EFTCAMBModelInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBModelFeedback(self, print_params)
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model) :: self         !< the base class
    logical, optional    :: print_params !< optional flag that decised whether to print numerical values
    !! of the parameters.
    end subroutine EFTCAMBModelFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model.
    subroutine EFTCAMBModelParameterNames(self, i, name)
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model) :: self   !< the base class
    integer, intent(in)  :: i      !< the index of the parameter
    character(*), intent(out) :: name   !< the output name of the i-th parameter
    end subroutine EFTCAMBModelParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model.
    subroutine EFTCAMBModelParameterNamesLatex(self, i, latexname)
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model) :: self       !< the base class
    integer, intent(in)  :: i          !< The index of the parameter
    character(*), intent(out) :: latexname  !< the output latex name of the i-th parameter
    end subroutine EFTCAMBModelParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model.
    subroutine EFTCAMBModelParameterValues(self, i, value)
    use precision
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model) :: self   !< the base class
    integer, intent(in)  :: i      !< The index of the parameter
    real(dl), intent(out) :: value  !< the output value of the i-th parameter
    end subroutine EFTCAMBModelParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBModelBackgroundEFTFunctions(self, a, eft_par_cache, eft_cache)
    use precision
    use EFTCAMB_cache
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model)                          :: self          !< the base class.
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    end subroutine EFTCAMBModelBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the second order EFT functions at a given time.
    subroutine EFTCAMBModelSecondOrderEFTFunctions(self, a, eft_par_cache, eft_cache)
    use precision
    use EFTCAMB_cache
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model)                          :: self          !< the base class.
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    end subroutine EFTCAMBModelSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).
    function EFTCAMBModelComputeDtauda(self, a, eft_par_cache, eft_cache)
    use precision
    use EFTCAMB_cache
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model)                          :: self          !< the base class.
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    real(dl)                                     :: EFTCAMBModelComputeDtauda !< the output dtauda
    end function EFTCAMBModelComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBModelComputeAdotoa(self, a, eft_par_cache, eft_cache)
    use precision
    use EFTCAMB_cache
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model)                          :: self          !< the base class.
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    end subroutine EFTCAMBModelComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBModelComputeHubbleDer(self, a, eft_par_cache, eft_cache)
    use precision
    use EFTCAMB_cache
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model)                          :: self          !< the base class.
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    end subroutine EFTCAMBModelComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes additional derivatives of the EFT functions, as needed by the
    !> stability calculation.
    subroutine EFTCAMBModelComputeAdditionalDer(self, a, eft_par_cache, eft_cache)
    use precision
    use EFTCAMB_cache
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model)                          :: self          !< the base class.
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    end subroutine EFTCAMBModelComputeAdditionalDer

    ! ---------------------------------------------------------------------------------------------

    end interface

    contains

    ! ---------------------------------------------------------------------------------------------
    ! EFTCAMB abstract model implementation: the following are all the procedures that can be
    ! be safely implemented for the abstract class and are not harmful if not overritten.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the name and latex name of the model.
    subroutine EFTCAMBModelInitialize(self, name, latexname)

    implicit none

    class(EFTCAMB_model)     :: self      !< the base class
    character(*), intent(in) :: name      !< the name of the function
    character(*), intent(in) :: latexname !< the latex name of the function

    self%name = trim(name)
    self%name_latex = trim(latexname)

    end subroutine EFTCAMBModelInitialize

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of the model, if needed.
    subroutine EFTCAMBModelInitBackground(self, params_cache, feedback_level, success, outroot)

    implicit none

    class(EFTCAMB_model)                          :: self           !< the base class
    !< a EFTCAMB parameter cache containing cosmological parameters
    type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache
    integer, intent(in)  :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
    logical, intent(out) :: success        !< wether the background initialization succeded or not
    character(len=*), intent(in), optional :: outroot        !< the output root for the debug files

    success = .true.

    end subroutine EFTCAMBModelInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes \rho_Q and P_Q. For details refer to the numerical notes.
    subroutine EFTCAMBModelComputeRhoQPQ(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: a2, adotoa2, aomegaP

    ! precompute some parts:
    a2 = a*a
    adotoa2 = eft_cache%adotoa**2
    aomegaP = a*eft_cache%EFTOmegaP

    ! do the computations:
    eft_cache%grhoq = 2._dl*eft_cache%EFTc - eft_cache%EFTLambda - 3._dl*adotoa2*aomegaP
    eft_cache%gpresq = eft_cache%EFTLambda + a2*adotoa2*eft_cache%EFTOmegaPP + aomegaP*(eft_cache%Hdot + 2._dl*adotoa2)
    eft_cache%grhodotq = 3._dl*eft_cache%adotoa*(-eft_cache%grhoq - eft_cache%gpresq + adotoa2*aomegaP)
    eft_cache%gpresdotq = eft_cache%EFTLambdadot &
        + adotoa2*eft_cache%adotoa*(a*a2*eft_cache%EFTOmegaPPP - 2._dl*aomegaP + 2._dl*a2*eft_cache%EFTOmegaPP) &
        + aomegaP*eft_cache%Hdotdot &
        + 3._dl*eft_cache%adotoa*eft_cache%Hdot*(aomegaP + a2*eft_cache%EFTOmegaPP)

    end subroutine EFTCAMBModelComputeRhoQPQ

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the Einstein equations factors. For details refer to the numerical notes.
    !! The implementation might look a bit intricated because this is a crucial part performance wise.
    subroutine EFTCAMBModelComputeEinsteinFactors(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: one_plus_omega, adotoa2, aomegaP, k2, a2, hdot_m_adotoa2, pidot_p_H_pi

    ! precompute some common parts:
    one_plus_omega = 1._dl + eft_cache%EFTOmegaV
    adotoa2 = eft_cache%adotoa**2
    hdot_m_adotoa2 = eft_cache%Hdot - adotoa2
    aomegaP = a*eft_cache%EFTOmegaP
    k2 = eft_cache%k**2
    a2 = a**2
    pidot_p_H_pi = eft_cache%pidot + eft_cache%adotoa*eft_cache%pi

    ! compute the coefficients:
    eft_cache%EFTeomF = 1.5_dl/(eft_cache%k*one_plus_omega)* &
        ((eft_cache%grhoq + eft_cache%gpresq)*eft_cache%pi + &
        (aomegaP*eft_cache%adotoa + a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)*pidot_p_H_pi + &
        eft_cache%pi*(k2*(eft_cache%EFTGamma3V + eft_cache%EFTGamma4V) - &
        (3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V)*hdot_m_adotoa2))
    eft_cache%EFTeomG = +1._dl + &
        0.5/one_plus_omega*(aomegaP + a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V/eft_cache%adotoa + &
        3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V)
    eft_cache%EFTeomL = +0.5_dl/one_plus_omega*((2._dl*eft_cache%EFTc*pidot_p_H_pi + eft_cache%grhodotq*eft_cache%pi)/eft_cache%adotoa &
        - 3._dl*aomegaP*((3._dl*adotoa2 - eft_cache%Hdot + k2/3._dl)*eft_cache%pi + eft_cache%adotoa*eft_cache%pidot) &
        + 4._dl*a2*eft_par_cache%h0_Mpc**2*eft_cache%EFTGamma1V/eft_cache%adotoa*pidot_p_H_pi &
        + 3._dl*a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V*((eft_cache%Hdot - 2._dl*adotoa2 - k2/3._dl)*eft_cache%pi/eft_cache%adotoa - eft_cache%pidot) &                                                                                  ! Gamma2
        - eft_cache%pi*(k2 - 3._dl*hdot_m_adotoa2)*(3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V) + 8._dl*eft_cache%EFTGamma6V*k2*pidot_p_H_pi/eft_cache%adotoa) &
        + 2._dl*eft_cache%EFTGamma5V/one_plus_omega*eft_cache%k**2*eft_cache%pi
    eft_cache%EFTeomM = eft_cache%gpresdotq*eft_cache%pi + &
        (eft_cache%grhoq + eft_cache%gpresq + a2*adotoa2*eft_cache%EFTOmegaPP)*pidot_p_H_pi + &
        aomegaP*eft_cache%adotoa*(eft_cache%pidotdot + &
        (eft_cache%Hdot + 4._dl*adotoa2)*eft_cache%pidot/eft_cache%adotoa + &
        2._dl*(eft_cache%Hdot + 3._dl*adotoa2 + k2/3._dl)*eft_cache%pi) + &
        a*eft_par_cache%h0_Mpc*(eft_cache%EFTGamma2V*eft_cache%pidotdot + &
        (4._dl*eft_cache%EFTGamma2V + a*eft_cache%EFTGamma2P)*eft_cache%adotoa*eft_cache%pidot + &
        (3._dl*adotoa2*eft_cache%EFTGamma2V + eft_cache%Hdot*eft_cache%EFTGamma2V + a*adotoa2*eft_cache%EFTGamma2P)* &
        eft_cache%pi) - (hdot_m_adotoa2 - k2/3._dl)* &
        ((3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V)*eft_cache%pidot + &
        2._dl*eft_cache%adotoa*(+3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + 1.5_dl*a*eft_cache%EFTGamma3P + &
        0.5_dl*a*eft_cache%EFTGamma4P)*eft_cache%pi) - &
        (3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V)* &
        (eft_cache%Hdotdot - 2._dl*eft_cache%adotoa*eft_cache%Hdot)*eft_cache%pi - &
        4._dl*eft_cache%EFTGamma5V*k2*pidot_p_H_pi/3._dl
    eft_cache%EFTeomN = eft_cache%k/one_plus_omega* &
        (eft_cache%adotoa*eft_cache%pi*(-aomegaP + 2._dl*eft_cache%EFTGamma4V + a*eft_cache%EFTGamma4P) + &
        eft_cache%EFTGamma4V*eft_cache%pidot + 2._dl*eft_cache%EFTGamma5V*pidot_p_H_pi)
    eft_cache%EFTeomNdot = eft_cache%k/one_plus_omega* &
        (-eft_cache%Hdot*aomegaP*eft_cache%pi - eft_cache%adotoa*aomegaP*eft_cache%pidot - &
        adotoa2*(aomegaP + a2*eft_cache%EFTOmegaPP - aomegaP**2/one_plus_omega)*eft_cache%pi + &
        eft_cache%EFTGamma4V*eft_cache%pidotdot + &
        a*eft_cache%adotoa*eft_cache%pidot* &
        (+eft_cache%EFTGamma4P - eft_cache%EFTGamma4V*eft_cache%EFTOmegaP/one_plus_omega) + &
        2._dl*(eft_cache%EFTGamma4V + 0.5_dl*a*eft_cache%EFTGamma4P)* &
        (eft_cache%Hdot*eft_cache%pi + eft_cache%adotoa*eft_cache%pidot) + &
        2._dl*a*adotoa2*eft_cache%pi* &
        (+0.5_dl*a*eft_cache%EFTGamma4PP + 1.5_dl*eft_cache%EFTGamma4P - &
        eft_cache%EFTOmegaP/one_plus_omega*(eft_cache%EFTGamma4V + 0.5_dl*a*eft_cache%EFTGamma4P)) + &
        2._dl*eft_cache%EFTGamma5V* &
        (eft_cache%pidotdot + eft_cache%adotoa*eft_cache%pidot + eft_cache%Hdot*eft_cache%pi) + &
        2._dl*eft_cache%adotoa*pidot_p_H_pi*(+a*eft_cache%EFTGamma5P - eft_cache%EFTGamma5V*aomegaP/one_plus_omega))
    eft_cache%EFTeomU = 1._dl + (+1.5_dl*eft_cache%EFTGamma3V + 0.5_dl*eft_cache%EFTGamma4V)/one_plus_omega
    eft_cache%EFTeomV = +0.5_dl/one_plus_omega*(aomegaP - 2._dl*eft_cache%EFTGamma4V - a*eft_cache%EFTGamma4P)
    eft_cache%EFTeomVdot = 0.5_dl*eft_cache%adotoa/one_plus_omega* &
        (aomegaP - 3._dl*a*eft_cache%EFTGamma4P + a2*(eft_cache%EFTOmegaPP - eft_cache%EFTGamma4PP) + &
        aomegaP/one_plus_omega*(-aomegaP + 2._dl*eft_cache%EFTGamma4V + a*eft_cache%EFTGamma4P))
    eft_cache%EFTeomX = 1._dl - eft_cache%EFTGamma4V/one_plus_omega
    eft_cache%EFTeomXdot = -a*eft_cache%adotoa/one_plus_omega*(+eft_cache%EFTGamma4P &
        - eft_cache%EFTGamma4V*eft_cache%EFTOmegaP/one_plus_omega)
    eft_cache%EFTeomY = +0.5_dl/one_plus_omega*(aomegaP &
        + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V &
        + 0.5_dl*a*(3._dl*eft_cache%EFTGamma3P + eft_cache%EFTGamma4P))
    eft_cache%EFTeomQ = 1._dl + 2._dl*eft_cache%EFTGamma5V/one_plus_omega

    end subroutine EFTCAMBModelComputeEinsteinFactors

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the pi field equations factors. For details refer to the numerical notes.
    subroutine EFTCAMBModelComputePiFactors(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: one_plus_omega, adotoa2, aomegaP, k2, a2, hdot_m_adotoa2, adotoa02

    ! precompute some common parts:
    k2 = eft_cache%k**2
    a2 = a**2
    one_plus_omega = 1._dl + eft_cache%EFTOmegaV
    adotoa2 = eft_cache%adotoa**2
    adotoa02 = eft_par_cache%h0_Mpc**2
    hdot_m_adotoa2 = eft_cache%Hdot - adotoa2
    aomegaP = a*eft_cache%EFTOmegaP

    ! compute the coefficients:
    eft_cache%EFTpiA1 = eft_cache%EFTc + 2._dl*a2*adotoa02*eft_cache%EFTGamma1V + &
        1.5_dl*a2*(eft_cache%adotoa*eft_cache%EFTOmegaP + eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)**2/ &
        (2._dl*one_plus_omega + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V)
    !
    eft_cache%EFTpiA2 = +4._dl*eft_cache%EFTGamma6V
    !
    eft_cache%EFTpiB1 = eft_cache%EFTcdot + 4._dl*eft_cache%adotoa*eft_cache%EFTc + &
        8._dl*a2*eft_cache%adotoa*adotoa02*(eft_cache%EFTGamma1V + 0.25_dl*a*eft_cache%EFTGamma1P) - &
        a*(eft_cache%adotoa*eft_cache%EFTOmegaP + eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)/ &
        (4._dl*one_plus_omega + 6._dl*eft_cache%EFTGamma3V + 2._dl*eft_cache%EFTGamma4V)* &
        (-3._dl*(eft_cache%grhoq + eft_cache%gpresq) - 3._dl*adotoa2*aomegaP*(4._dl + eft_cache%Hdot/adotoa2) - &
        3._dl*a2*adotoa2*eft_cache%EFTOmegaPP - &
        3._dl*a*eft_cache%adotoa*eft_par_cache%h0_Mpc*(4._dl*eft_cache%EFTGamma2V + a*eft_cache%EFTGamma2P) + &
        (9._dl*eft_cache%EFTGamma3V + 3._dl*eft_cache%EFTGamma4V)*hdot_m_adotoa2) + &
        1._dl/(one_plus_omega + 2._dl*eft_cache%EFTGamma5V)* &
        (eft_cache%adotoa*aomegaP + 2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V + a*eft_cache%EFTGamma5P) - &
        one_plus_omega*(eft_cache%adotoa*aomegaP + a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)/ &
        (2._dl*one_plus_omega + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))* &
        (-eft_cache%EFTc + 1.5_dl*a*adotoa2*eft_cache%EFTOmegaP - &
        2._dl*a2*eft_par_cache%h0_Mpc**2*eft_cache%EFTGamma1V + &
        1.5_dl*a*eft_cache%adotoa*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)
    !
    eft_cache%EFTpiB2 = +4._dl*eft_cache%adotoa*(2._dl*eft_cache%EFTGamma6V + a*eft_cache%EFTGamma6P) + &
        a*(eft_cache%EFTGamma4V + 2._dl*eft_cache%EFTGamma5V)/ &
        (2._dl*(1._dl + eft_cache%EFTOmegaV) - 2._dl*eft_cache%EFTGamma4V)* &
        (eft_cache%adotoa*eft_cache%EFTOmegaP + eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V) - &
        a*(eft_cache%adotoa*eft_cache%EFTOmegaP + eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)/ &
        (4._dl*one_plus_omega + 6._dl*eft_cache%EFTGamma3V + 2._dl*eft_cache%EFTGamma4V)* &
        (+(-3._dl*eft_cache%EFTGamma3V - eft_cache%EFTGamma4V + 4._dl*eft_cache%EFTGamma5V)) + &
        1._dl/(1._dl + eft_cache%EFTOmegaV + 2._dl*eft_cache%EFTGamma5V)* &
        (a*eft_cache%adotoa*eft_cache%EFTOmegaP + &
        2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V + a*eft_cache%EFTGamma5P) - &
        (1._dl + eft_cache%EFTOmegaV)* &
        (a*eft_cache%adotoa*eft_cache%EFTOmegaP + a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)/ &
        (2._dl*one_plus_omega + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))*(-4._dl*eft_cache%EFTGamma6V)
    !
    eft_cache%EFTpiC = +eft_cache%adotoa*eft_cache%EFTcdot + (6._dl*adotoa2 - 2._dl*eft_cache%Hdot)*eft_cache%EFTc + &
        1.5_dl*a*eft_cache%adotoa*eft_cache%EFTOmegaP*(eft_cache%Hdotdot - 2._dl*eft_cache%adotoa**3) + &
        6._dl*(a*eft_cache%adotoa*eft_par_cache%h0_Mpc)**2*eft_cache%EFTGamma1V + &
        2._dl*a2*eft_cache%Hdot*adotoa02*eft_cache%EFTGamma1V + 2._dl*a**3*adotoa2*adotoa02*eft_cache%EFTGamma1P + &
        1.5_dl*hdot_m_adotoa2**2*(eft_cache%EFTGamma4V + 3._dl*eft_cache%EFTGamma3V) + &
        4.5_dl*eft_cache%adotoa*eft_par_cache%h0_Mpc*a*hdot_m_adotoa2* &
        (eft_cache%EFTGamma2V + a*eft_cache%EFTGamma2P/3._dl) + &
        0.5_dl*a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V* &
        (3._dl*eft_cache%Hdotdot - 12._dl*eft_cache%Hdot*eft_cache%adotoa + 6._dl*eft_cache%adotoa**3) - &
        a*(eft_cache%adotoa*eft_cache%EFTOmegaP + eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)/ &
        (4._dl*(1._dl + eft_cache%EFTOmegaV) + 6._dl*eft_cache%EFTGamma3V + 2._dl*eft_cache%EFTGamma4V)* &
        (-3._dl*eft_cache%gpresdotq - 3._dl*eft_cache%adotoa*(eft_cache%grhoq + eft_cache%gpresq) - &
        3._dl*a*eft_cache%adotoa**3*(a*eft_cache%EFTOmegaPP + 6._dl*eft_cache%EFTOmegaP) - &
        6._dl*a*eft_cache%adotoa*eft_cache%Hdot*eft_cache%EFTOmegaP + &
        3._dl*(eft_cache%Hdotdot - 2._dl*eft_cache%adotoa*eft_cache%Hdot)* &
        (eft_cache%EFTGamma4V + 3._dl*eft_cache%EFTGamma3V) + &
        6._dl*eft_cache%adotoa*hdot_m_adotoa2* &
        (3._dl*eft_cache%EFTGamma3V + 1.5_dl*a*eft_cache%EFTGamma3P + eft_cache%EFTGamma4V + &
        0.5_dl*a*eft_cache%EFTGamma4P) - &
        3._dl*a*eft_par_cache%h0_Mpc* &
        (3._dl*adotoa2*eft_cache%EFTGamma2V + eft_cache%Hdot*eft_cache%EFTGamma2V + a*adotoa2*eft_cache%EFTGamma2P)) + &
        1._dl/(1._dl + eft_cache%EFTOmegaV + 2._dl*eft_cache%EFTGamma5V)* &
        (a*eft_cache%adotoa*eft_cache%EFTOmegaP + &
        2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V + a*eft_cache%EFTGamma5P) - &
        (1._dl + eft_cache%EFTOmegaV)* &
        (a*eft_cache%adotoa*eft_cache%EFTOmegaP + a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)/ &
        (2._dl*(1._dl + eft_cache%EFTOmegaV) + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))* &
        (-0.5*eft_cache%grhodotq - eft_cache%adotoa*eft_cache%EFTc + &
        1.5_dl*a*eft_cache%adotoa*eft_cache%EFTOmegaP*(3._dl*adotoa2 - eft_cache%Hdot) - &
        2._dl*a2*eft_cache%adotoa*adotoa02*eft_cache%EFTGamma1V - &
        1.5_dl*a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V*(eft_cache%Hdot - 2._dl*adotoa2) - &
        3._dl*eft_cache%adotoa*hdot_m_adotoa2*(1.5_dl*eft_cache%EFTGamma3V + 0.5_dl*eft_cache%EFTGamma4V))
    !
    eft_cache%EFTpiD1 = eft_cache%EFTc - &
        0.5_dl*a*eft_cache%adotoa*eft_par_cache%h0_Mpc*(eft_cache%EFTGamma2V + a*eft_cache%EFTGamma2P) - &
        hdot_m_adotoa2*(3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V) + &
        4._dl*(eft_cache%Hdot*eft_cache%EFTGamma6V + adotoa2*eft_cache%EFTGamma6V + a*adotoa2*eft_cache%EFTGamma6P) + &
        2._dl*(eft_cache%Hdot*eft_cache%EFTGamma5V + a*adotoa2*eft_cache%EFTGamma5P) - &
        a*(eft_cache%adotoa*eft_cache%EFTOmegaP + eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)/ &
        (4._dl*(1._dl + eft_cache%EFTOmegaV) + 6._dl*eft_cache%EFTGamma3V + 2._dl*eft_cache%EFTGamma4V)* &
        (-2._dl*a*eft_cache%adotoa*eft_cache%EFTOmegaP + 4._dl*eft_cache%adotoa*eft_cache%EFTGamma5V - &
        2._dl*eft_cache%adotoa*(3._dl*eft_cache%EFTGamma3V + 1.5_dl*a*eft_cache%EFTGamma3P + eft_cache%EFTGamma4V + &
        0.5_dl*a*eft_cache%EFTGamma4P)) + &
        1._dl/(1._dl + eft_cache%EFTOmegaV + 2._dl*eft_cache%EFTGamma5V)* &
        (a*eft_cache%adotoa*eft_cache%EFTOmegaP + &
        2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V + a*eft_cache%EFTGamma5P) - &
        (1._dl + eft_cache%EFTOmegaV)* &
        (a*eft_cache%adotoa*eft_cache%EFTOmegaP + a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)/ &
        (2._dl*(1._dl + eft_cache%EFTOmegaV) + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))* &
        (+0.5_dl*a*eft_cache%adotoa*eft_cache%EFTOmegaP - 2._dl*eft_cache%adotoa*eft_cache%EFTGamma5V + &
        0.5_dl*a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V + 1.5_dl*eft_cache%adotoa*eft_cache%EFTGamma3V + &
        0.5_dl*eft_cache%adotoa*eft_cache%EFTGamma4V - 4._dl*eft_cache%adotoa*eft_cache%EFTGamma6V) + &
        (eft_cache%EFTGamma4V + 2._dl*eft_cache%EFTGamma5V)/ &
        (2._dl*(1._dl + eft_cache%EFTOmegaV) - 2._dl*eft_cache%EFTGamma4V)* &
        (eft_cache%grhoq + eft_cache%gpresq + a*adotoa2*eft_cache%EFTOmegaP - eft_cache%EFTGamma4V*hdot_m_adotoa2 + &
        a*eft_cache%adotoa*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V - 3._dl*eft_cache%EFTGamma3V*hdot_m_adotoa2)
    !
    eft_cache%EFTpiD2 = +(+0.5_dl*eft_cache%EFTGamma3V + 0.5_dl*eft_cache%EFTGamma4V + &
        (eft_cache%EFTGamma4V + 2._dl*eft_cache%EFTGamma5V)/ &
        (2._dl*(1._dl + eft_cache%EFTOmegaV) - 2._dl*eft_cache%EFTGamma4V)* &
        (eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))
    !
    eft_cache%EFTpiE = (eft_cache%EFTc - 1.5_dl*a*adotoa2*eft_cache%EFTOmegaP - &
        0.5_dl*a*eft_cache%adotoa*eft_par_cache%h0_Mpc*(2._dl*eft_cache%EFTGamma2V + a*eft_cache%EFTGamma2P) + &
        0.5_dl*eft_cache%EFTGamma3V*(k2 - 3._dl*eft_cache%Hdot + 3._dl*adotoa2) + &
        0.5_dl*eft_cache%EFTGamma4V*(k2 - eft_cache%Hdot + adotoa2) - &
        a*(eft_cache%adotoa*eft_cache%EFTOmegaP + eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)/ &
        (4._dl*(1._dl + eft_cache%EFTOmegaV) + 6._dl*eft_cache%EFTGamma3V + 2._dl*eft_cache%EFTGamma4V)* &
        (-2._dl*eft_cache%adotoa*(a*eft_cache%EFTOmegaP + 2._dl*(1._dl + eft_cache%EFTOmegaV)) - &
        2._dl*eft_cache%adotoa*(3._dl*eft_cache%EFTGamma3V + 1.5_dl*a*eft_cache%EFTGamma3P + eft_cache%EFTGamma4V + &
        0.5_dl*a*eft_cache%EFTGamma4P)) + &
        1._dl/(1._dl + eft_cache%EFTOmegaV + 2._dl*eft_cache%EFTGamma5V)* &
        (a*eft_cache%adotoa*eft_cache%EFTOmegaP + &
        2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V + a*eft_cache%EFTGamma5P) - &
        (1._dl + eft_cache%EFTOmegaV)* &
        (a*eft_cache%adotoa*eft_cache%EFTOmegaP + a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)/ &
        (2._dl*(1._dl + eft_cache%EFTOmegaV) + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))* &
        (+eft_cache%adotoa*(1._dl + eft_cache%EFTOmegaV + 0.5_dl*a*eft_cache%EFTOmegaP) + &
        0.5_dl*a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V + 1.5_dl*eft_cache%adotoa*eft_cache%EFTGamma3V + &
        0.5_dl*eft_cache%adotoa*eft_cache%EFTGamma4V) + &
        (eft_cache%EFTGamma4V + 2._dl*eft_cache%EFTGamma5V)/ &
        (2._dl*(1._dl + eft_cache%EFTOmegaV) - 2._dl*eft_cache%EFTGamma4V)*k2* &
        (eft_cache%EFTGamma4V + eft_cache%EFTGamma3V))*eft_cache%k*eft_cache%z + &
        1._dl*a*(eft_cache%adotoa*eft_cache%EFTOmegaP + eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)/ &
        (4._dl*(1._dl + eft_cache%EFTOmegaV) + 6._dl*eft_cache%EFTGamma3V + 2._dl*eft_cache%EFTGamma4V)* &
        (eft_cache%grhog_t*eft_cache%clxg + eft_cache%grhor_t*eft_cache%clxr + 3._dl*eft_cache%dgpnu) + &
        (eft_cache%EFTGamma4V + 2._dl*eft_cache%EFTGamma5V)/ &
        (2._dl*(1._dl + eft_cache%EFTOmegaV) - 2._dl*eft_cache%EFTGamma4V)*eft_cache%k*eft_cache%dgq - &
        0.5_dl/(1._dl + eft_cache%EFTOmegaV + 2._dl*eft_cache%EFTGamma5V)* &
        (a*eft_cache%adotoa*eft_cache%EFTOmegaP + &
        2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V + a*eft_cache%EFTGamma5P) - &
        (1._dl + eft_cache%EFTOmegaV)* &
        (a*eft_cache%adotoa*eft_cache%EFTOmegaP + a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V)/ &
        (2._dl*(1._dl + eft_cache%EFTOmegaV) + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))*eft_cache%dgrho

    end subroutine EFTCAMBModelComputePiFactors

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the factors for the tensor propagation equation. For details refer to the numerical notes.
    subroutine EFTCAMBModelComputeTensorFactors(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: one_plus_omega

    ! precompute some common parts:
    one_plus_omega = 1._dl + eft_cache%EFTOmegaV

    ! compute the coefficients:
    eft_cache%EFTAT = one_plus_omega - eft_cache%EFTGamma4V
    eft_cache%EFTBT = 2._dl*eft_cache%adotoa* &
        (one_plus_omega - eft_cache%EFTGamma4V + 0.5_dl*a*eft_cache%EFTOmegaP - 0.5_dl*a*eft_cache%EFTGamma4P)
    eft_cache%EFTDT = one_plus_omega

    end subroutine EFTCAMBModelComputeTensorFactors

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the kinetic and gradient terms. For details refer to the numerical notes.
    subroutine EFTCAMBModelComputeStabilityFactors(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: F1, F2, F3, F1dot, F2dot, rho_d, rho_d_dot, F1dotdot, F2dotdot, F3dot, F3dotdot
    real(dl) :: EFTcP, EFTcPP

    eft_cache%EFT_kinetic = 9._dl*(1._dl + eft_cache%EFTOmegaV - eft_cache%EFTGamma4V)* &
        (4._dl*eft_cache%EFTc*(1._dl + eft_cache%EFTOmegaV - eft_cache%EFTGamma4V) + &
        3._dl*eft_cache%adotoa**2*eft_cache%EFTOmegaP**2*a**2 + &
        a**2*eft_par_cache%h0_Mpc*(eft_par_cache%h0_Mpc* &
        (3._dl*eft_cache%EFTGamma2V**2 + &
        8._dl*eft_cache%EFTGamma1V*(1._dl + eft_cache%EFTOmegaV - eft_cache%EFTGamma4V)) + &
        6._dl*eft_cache%adotoa*eft_cache%EFTGamma2V*eft_cache%EFTOmegaP))
    !
    eft_cache%EFT_gradient = 9._dl* &
        (8._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5P - 16._dl*eft_cache%adotoa**2*eft_cache%EFTGamma5V**2 + &
        16._dl*eft_cache%EFTc*eft_cache%EFTGamma5V**2 - &
        2._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma4P*eft_cache%EFTOmegaP + &
        4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaP - &
        4._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP - &
        4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma4P*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP - &
        8._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaP + &
        3._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTOmegaP**2 + &
        4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP**2 + &
        16._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaV + &
        16._dl*eft_cache%EFTc*eft_cache%EFTGamma5V*eft_cache%EFTOmegaV - &
        16._dl*eft_cache%adotoa**2*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaV - &
        2._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma4P*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV + &
        4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV - &
        4._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV + &
        3._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTOmegaP**2*eft_cache%EFTOmegaV + &
        8._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaV**2 - &
        a**2*eft_cache%EFTGamma2V**2*eft_par_cache%h0_Mpc**2*(1 + eft_cache%EFTOmegaV) + &
        4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaPP* &
        (1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) + &
        4._dl*eft_cache%EFTc*(4._dl*eft_cache%EFTGamma5V + (1._dl + eft_cache%EFTOmegaV)**2) - &
        2._dl*a*eft_cache%adotoa*eft_par_cache%h0_Mpc* &
        (a*eft_cache%EFTGamma2P*(1._dl - eft_cache%EFTGamma4V + eft_cache%EFTOmegaV)* &
        (1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) + &
        eft_cache%EFTGamma2V*(eft_cache%EFTGamma4V* &
        (-1._dl + 2._dl*a*eft_cache%EFTGamma5P + 2._dl*eft_cache%EFTGamma5V + a*eft_cache%EFTOmegaP - &
        eft_cache%EFTOmegaV) + (1._dl + eft_cache%EFTOmegaV)* &
        (1._dl + a*(eft_cache%EFTGamma4P - 2._dl*eft_cache%EFTGamma5P - eft_cache%EFTOmegaP) + eft_cache%EFTOmegaV) - &
        2._dl*eft_cache%EFTGamma5V*(1._dl - a*eft_cache%EFTGamma4P + a*eft_cache%EFTOmegaP + eft_cache%EFTOmegaV))) + &
        8._dl*eft_cache%EFTGamma5V*eft_cache%Hdot + 16._dl*eft_cache%EFTGamma5V**2*eft_cache%Hdot + &
        4._dl*a*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP*eft_cache%Hdot + &
        8._dl*a*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaP*eft_cache%Hdot + &
        16._dl*eft_cache%EFTGamma5V*eft_cache%EFTOmegaV*eft_cache%Hdot + &
        16._dl*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaV*eft_cache%Hdot + &
        4._dl*a*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV*eft_cache%Hdot + &
        8._dl*eft_cache%EFTGamma5V*eft_cache%EFTOmegaV**2*eft_cache%Hdot + &
        4._dl*eft_cache%EFTGamma4V**2* &
        (eft_cache%adotoa**2*(1._dl + 2._dl*a*eft_cache%EFTGamma5P + 4._dl*eft_cache%EFTGamma5V + &
        a*eft_cache%EFTOmegaP + eft_cache%EFTOmegaV) - &
        (1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV)*eft_cache%Hdot) + &
        2._dl*eft_cache%EFTGamma4V* &
        (eft_cache%adotoa**2*(-(a**2*eft_cache%EFTOmegaP**2) + &
        a**2*eft_cache%EFTOmegaPP*(1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) - &
        4._dl*(1._dl + eft_cache%EFTOmegaV)* &
        (1._dl + 2._dl*a*eft_cache%EFTGamma5P + 4._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) - &
        a*eft_cache%EFTOmegaP*(3._dl + 2._dl*a*eft_cache%EFTGamma5P + 2._dl*eft_cache%EFTGamma5V + &
        3._dl*eft_cache%EFTOmegaV)) + &
        (1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV)* &
        (4._dl + a*eft_cache%EFTOmegaP + 4._dl*eft_cache%EFTOmegaV)*eft_cache%Hdot))

    end subroutine EFTCAMBModelComputeStabilityFactors

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the scales of the scalar no-ghost (EFT_kinetic) and gradient (EFT_gradient)
    !! coefficients of EFTCAMBModelComputeStabilityFactors: the sums of the magnitudes of their terms.
    !! They are the reference for the relative numerical floor EFTCAMB_stability_threshold of the stability check:
    !! a coefficient that cancels to roundoff of its terms has no definite sign.
    subroutine EFTCAMBStabilityScales(a, eft_par_cache, eft_cache, kinetic_scale, gradient_scale)

    implicit none

    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    real(dl), intent(out) :: kinetic_scale    !< sum of the magnitudes of the terms of EFT_kinetic
    real(dl), intent(out) :: gradient_scale   !< sum of the magnitudes of the terms of EFT_gradient

    real(dl) :: M

    M = 1._dl + eft_cache%EFTOmegaV - eft_cache%EFTGamma4V
    kinetic_scale = 9._dl*abs(M)*(abs(4._dl*eft_cache%EFTc*M) + &
        3._dl*eft_cache%adotoa**2*eft_cache%EFTOmegaP**2*a**2 + &
        a**2*eft_par_cache%h0_Mpc**2*(3._dl*eft_cache%EFTGamma2V**2 + abs(8._dl*eft_cache%EFTGamma1V*M)) + &
        abs(6._dl*a**2*eft_par_cache%h0_Mpc*eft_cache%adotoa*eft_cache%EFTGamma2V*eft_cache%EFTOmegaP))
    !
    gradient_scale = 9._dl*(abs(8._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5P) + &
        abs(16._dl*eft_cache%adotoa**2*eft_cache%EFTGamma5V**2) + &
        abs(16._dl*eft_cache%EFTc*eft_cache%EFTGamma5V**2) + &
        abs(2._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma4P*eft_cache%EFTOmegaP) + &
        abs(4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaP) + &
        abs(4._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP) + &
        abs(4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma4P*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP) + &
        abs(8._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaP) + &
        abs(3._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTOmegaP**2) + &
        abs(4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP**2) + &
        abs(16._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaV) + &
        abs(16._dl*eft_cache%EFTc*eft_cache%EFTGamma5V*eft_cache%EFTOmegaV) + &
        abs(16._dl*eft_cache%adotoa**2*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaV) + &
        abs(2._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma4P*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV) + &
        abs(4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV) + &
        abs(4._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV) + &
        abs(3._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTOmegaP**2*eft_cache%EFTOmegaV) + &
        abs(8._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaV**2) + &
        abs(a**2*eft_cache%EFTGamma2V**2*eft_par_cache%h0_Mpc**2*(1 + eft_cache%EFTOmegaV)) + &
        abs(4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaPP* (1._dl + 2._dl* &
        eft_cache%EFTGamma5V + eft_cache%EFTOmegaV)) + &
        abs(4._dl*eft_cache%EFTc*(4._dl*eft_cache%EFTGamma5V + (1._dl + eft_cache%EFTOmegaV)**2)) + &
        abs(2._dl*a*eft_cache%adotoa*eft_par_cache%h0_Mpc* (a*eft_cache%EFTGamma2P*(1._dl - eft_cache%EFTGamma4V &
        + eft_cache%EFTOmegaV)* (1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) + eft_cache%EFTGamma2V* &
        (eft_cache%EFTGamma4V* (-1._dl + 2._dl*a*eft_cache%EFTGamma5P + 2._dl*eft_cache%EFTGamma5V + a* &
        eft_cache%EFTOmegaP - eft_cache%EFTOmegaV) + (1._dl + eft_cache%EFTOmegaV)* (1._dl + a* &
        (eft_cache%EFTGamma4P - 2._dl*eft_cache%EFTGamma5P - eft_cache%EFTOmegaP) + eft_cache%EFTOmegaV) - 2._dl* &
        eft_cache%EFTGamma5V*(1._dl - a*eft_cache%EFTGamma4P + a*eft_cache%EFTOmegaP + eft_cache%EFTOmegaV)))) + &
        abs(8._dl*eft_cache%EFTGamma5V*eft_cache%Hdot) + abs(16._dl*eft_cache%EFTGamma5V**2*eft_cache%Hdot) + &
        abs(4._dl*a*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP*eft_cache%Hdot) + &
        abs(8._dl*a*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaP*eft_cache%Hdot) + &
        abs(16._dl*eft_cache%EFTGamma5V*eft_cache%EFTOmegaV*eft_cache%Hdot) + &
        abs(16._dl*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaV*eft_cache%Hdot) + &
        abs(4._dl*a*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV*eft_cache%Hdot) + &
        abs(8._dl*eft_cache%EFTGamma5V*eft_cache%EFTOmegaV**2*eft_cache%Hdot) + &
        abs(4._dl*eft_cache%EFTGamma4V**2* (eft_cache%adotoa**2*(1._dl + 2._dl*a*eft_cache%EFTGamma5P + 4._dl* &
        eft_cache%EFTGamma5V + a*eft_cache%EFTOmegaP + eft_cache%EFTOmegaV) - (1._dl + 2._dl* &
        eft_cache%EFTGamma5V + eft_cache%EFTOmegaV)*eft_cache%Hdot)) + &
        abs(2._dl*eft_cache%EFTGamma4V* (eft_cache%adotoa**2*(-(a**2*eft_cache%EFTOmegaP**2) + a**2* &
        eft_cache%EFTOmegaPP*(1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) - 4._dl*(1._dl &
        + eft_cache%EFTOmegaV)* (1._dl + 2._dl*a*eft_cache%EFTGamma5P + 4._dl*eft_cache%EFTGamma5V &
        + eft_cache%EFTOmegaV) - a*eft_cache%EFTOmegaP*(3._dl + 2._dl*a*eft_cache%EFTGamma5P + 2._dl* &
        eft_cache%EFTGamma5V + 3._dl*eft_cache%EFTOmegaV)) + (1._dl + 2._dl*eft_cache%EFTGamma5V &
        + eft_cache%EFTOmegaV)* (4._dl + a*eft_cache%EFTOmegaP + 4._dl*eft_cache%EFTOmegaV)*eft_cache%Hdot)))

    end subroutine EFTCAMBStabilityScales

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the kinetic and gradient terms. For details refer to the numerical notes.
    subroutine EFTCAMBModelComputeMassFactors(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: F1, F2, F3, F1dot, F2dot, rho_d, rho_d_dot, F1dotdot, F2dotdot, F3dot, F3dotdot
    real(dl) :: EFTcP, EFTcPP

    !< Define useful quantities used for the mass stability
    rho_d = (eft_cache%grhob_t + eft_cache%grhoc_t)/a**2._dl
    rho_d_dot = -3._dl*eft_cache%adotoa/a*rho_d
    EFTcP = eft_cache%EFTcdot/(a*eft_cache%adotoa) + 2._dl*eft_cache%EFTc/a

    !< Define some EFT functions needed by stability
    EFTcPP = (eft_cache%EFTcdotdot/eft_cache%adotoa + (3._dl - eft_cache%Hdot/eft_cache%adotoa)*eft_cache%EFTcdot + &
        2._dl*eft_cache%adotoa*eft_cache%EFTc)/(a**2*eft_cache%adotoa**2)

    !< general expressions for Fs
    F1 = 2._dl*(1._dl + eft_cache%EFTOmegaV) + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V

    F2 = eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V + &
        (eft_cache%adotoa*(2._dl*(1._dl + eft_cache%EFTOmegaV) + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))/ &
        a + eft_cache%adotoa*eft_cache%EFTOmegaP

    F3 = (2._dl*eft_cache%EFTc)/a**2. + 4._dl*eft_par_cache%h0_Mpc**2._dl*eft_cache%EFTGamma1V - &
        (6._dl*eft_par_cache%h0_Mpc*eft_cache%adotoa*eft_cache%EFTGamma2V)/a - &
        (3._dl*eft_cache%adotoa**2.* &
        (2._dl*(1._dl + eft_cache%EFTOmegaV) + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))/a**2. - &
        (6._dl*eft_cache%adotoa**2.*eft_cache%EFTOmegaP)/a

    F1dot = eft_cache%adotoa*(2._dl*eft_cache%EFTOmegaP + 3._dl*eft_cache%EFTGamma3P + eft_cache%EFTGamma4P)

    F2dot = ((2._dl + 2._dl*eft_cache%EFTOmegaV + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + &
        a*eft_cache%EFTOmegaP)*eft_cache%Hdot + &
        a**2._dl*eft_par_cache%h0_Mpc**2._dl*eft_cache%adotoa*eft_cache%EFTGamma2P + &
        eft_cache%adotoa**2._dl*(-2._dl - 2._dl*eft_cache%EFTOmegaV - 3._dl*eft_cache%EFTGamma3V - &
        eft_cache%EFTGamma4V + 2._dl*a*eft_cache%EFTOmegaP + 3._dl*a*eft_cache%EFTGamma3P + a*eft_cache%EFTGamma4P + &
        a**2._dl*eft_cache%EFTOmegaPP))/a**2._dl

    F3dot = -((4._dl*eft_cache%EFTc*eft_cache%adotoa + &
        6._dl*a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V*eft_cache%Hdot - &
        2._dl*eft_cache%adotoa*(a*EFTcP - &
        3._dl*(2._dl + 2._dl*eft_cache%EFTOmegaV + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + &
        2._dl*a*eft_cache%EFTOmegaP)*eft_cache%Hdot + 2._dl*a**3._dl*eft_par_cache%h0_Mpc**2*eft_cache%EFTGamma1P) + &
        6._dl*a*eft_par_cache%h0_Mpc*eft_cache%adotoa**2._dl*(-eft_cache%EFTGamma2V + a*eft_cache%EFTGamma2P) - &
        3._dl*eft_cache%adotoa**3._dl* &
        (4._dl + 4._dl*eft_cache%EFTOmegaV + 6._dl*eft_cache%EFTGamma3V + 2._dl*eft_cache%EFTGamma4V - &
        3._dl*a*eft_cache%EFTGamma3P - a*eft_cache%EFTGamma4P - 2._dl*a**2*eft_cache%EFTOmegaPP))/a**3._dl)

    F1dotdot = (2._dl*eft_cache%EFTOmegaP*eft_cache%Hdot + &
        eft_cache%Hdot*(3._dl*eft_cache%EFTGamma3P + eft_cache%EFTGamma4P) + &
        a*eft_cache%adotoa**2._dl*(2._dl*eft_cache%EFTOmegaPP + 3._dl*eft_cache%EFTGamma3PP + eft_cache%EFTGamma4PP))/a

    F2dotdot = (a**2._dl*eft_par_cache%h0_Mpc**2._dl*eft_cache%Hdot*eft_cache%EFTGamma2P + &
        eft_cache%adotoa*eft_cache%Hdot* &
        (-8._dl - 8._dl*eft_cache%EFTOmegaV - 12._dl*eft_cache%EFTGamma3V - 4._dl*eft_cache%EFTGamma4V + &
        5._dl*a*eft_cache%EFTOmegaP + 9._dl*a*eft_cache%EFTGamma3P + 3._dl*a*eft_cache%EFTGamma4P + &
        3._dl*a**2._dl*eft_cache%EFTOmegaPP) + &
        (2._dl + 2._dl*eft_cache%EFTOmegaV + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + &
        a*eft_cache%EFTOmegaP)*eft_cache%Hdotdot + &
        a**3._dl*eft_par_cache%h0_Mpc**2._dl*eft_cache%adotoa**2._dl*eft_cache%EFTGamma2PP + &
        eft_cache%adotoa**3._dl*(4._dl + 4._dl*eft_cache%EFTOmegaV + 6._dl*eft_cache%EFTGamma3V + &
        2._dl*eft_cache%EFTGamma4V - 4._dl*a*eft_cache%EFTOmegaP - 6._dl*a*eft_cache%EFTGamma3P - &
        2._dl*a*eft_cache%EFTGamma4P + 2._dl*a**2._dl*eft_cache%EFTOmegaPP + 3._dl*a**2._dl*eft_cache%EFTGamma3PP + &
        a**2._dl*eft_cache%EFTGamma4PP + a**3._dl*eft_cache%EFTOmegaPPP))/a**3._dl

    F3dotdot = (4._dl*eft_cache%EFTc*(3._dl*eft_cache%adotoa**2._dl - eft_cache%Hdot) + &
        2._dl*(a*EFTcP*eft_cache%Hdot - &
        3._dl*(2._dl + 2._dl*eft_cache%EFTOmegaV + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + &
        2*a*eft_cache%EFTOmegaP)*eft_cache%Hdot**2 + &
        2._dl*a**3._dl*eft_par_cache%h0_Mpc**2._dl*eft_cache%Hdot*eft_cache%EFTGamma1P - &
        3._dl*a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V*eft_cache%Hdotdot) - &
        6._dl*eft_cache%adotoa*(a*eft_par_cache%h0_Mpc*eft_cache%Hdot* &
        (-4._dl*eft_cache%EFTGamma2V + 3._dl*a*eft_cache%EFTGamma2P) + &
        (2._dl + 2._dl*eft_cache%EFTOmegaV + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + &
        2._dl*a*eft_cache%EFTOmegaP)*eft_cache%Hdotdot) + &
        eft_cache%adotoa**2._dl*(-8._dl*a*EFTcP + &
        3._dl*eft_cache%Hdot*(24._dl + 24._dl*eft_cache%EFTOmegaV + 36._dl*eft_cache%EFTGamma3V + &
        12._dl*eft_cache%EFTGamma4V + 4._dl*a*eft_cache%EFTOmegaP - 15._dl*a*eft_cache%EFTGamma3P - &
        5._dl*a*eft_cache%EFTGamma4P - 10._dl*a**2._dl*eft_cache%EFTOmegaPP) + &
        2._dl*a**2._dl*(EFTcPP + 2._dl*a**2._dl*eft_par_cache%h0_Mpc**2._dl*eft_cache%EFTGamma1PP)) - &
        6._dl*a*eft_par_cache%h0_Mpc*eft_cache%adotoa**3* &
        (2*eft_cache%EFTGamma2V + a*(-2._dl*eft_cache%EFTGamma2P + a*eft_cache%EFTGamma2PP)) - &
        3._dl*eft_cache%adotoa**4._dl* &
        (12._dl + 12._dl*eft_cache%EFTOmegaV + 18._dl*eft_cache%EFTGamma3V + 6._dl*eft_cache%EFTGamma4V - &
        4._dl*a*eft_cache%EFTOmegaP - 12._dl*a*eft_cache%EFTGamma3P - 4._dl*a*eft_cache%EFTGamma4P - &
        2._dl*a**2._dl*eft_cache%EFTOmegaPP + 3._dl*a**2._dl*eft_cache%EFTGamma3PP + a**2._dl*eft_cache%EFTGamma4PP + &
        2._dl*a**3._dl*eft_cache%EFTOmegaPPP))/a**4._dl

    !< compute the mass eigenvalues
    eft_cache%EFT_mu1 = ((-4._dl*F1*F2**2._dl*F3*(-2._dl*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2._dl)/ &
        (3._dl*F2**2._dl + F1*F3) + &
        (4._dl*F1*F3*(-2._dl*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2._dl)/(3._dl + (F1*F3)/F2**2._dl) - &
        (9._dl*F2**6*(1._dl + 1._dl/(3.*sqrt(F2**4._dl/(3._dl*F2**2._dl + 2._dl*F1*F3)**2._dl)))* &
        (6._dl*F1*F2*F3*(3._dl*F2**2._dl + F1*F3)*eft_cache%adotoa* &
        (F2*F3*F1dot + F1*(F2*F3dot + 2._dl*F3*(-F2dot + rho_d))) + &
        a*(-3._dl*F2**4._dl*F3**2._dl*F1dot**2._dl + &
        2._dl*F1*F2**2._dl*F3**2._dl*(-(F3*F1dot**2._dl) + 3._dl*F2**2._dl*F1dotdot) + &
        2._dl*F1**3._dl*F3*(-2._dl*F2**2._dl*F3dot**2._dl + F2*F3*(F2*F3dotdot - 2._dl*F3dot*(-2._dl*F2dot + rho_d)) - &
        2._dl*F3**2._dl*(F2dot**2._dl - F2dot*rho_d + rho_d**2._dl + F2*(F2dotdot - rho_d_dot))) + &
        F1**2._dl*F2*(4._dl*F3**3._dl*F1dot*F2dot + F2**3._dl*(-9._dl*F3dot**2._dl + 6._dl*F3*F3dotdot) + &
        2._dl*F2*F3**2._dl*(-(F1dot*F3dot) + F3*F1dotdot + 6._dl*(F2dot - rho_d)*rho_d) + &
        12._dl*F2**2._dl*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))))))/ &
        (a*(3._dl*F2**2._dl + F1*F3)**2._dl*(3._dl*F2**2._dl + 2._dl*F1*F3)) + &
        (2._dl*F1*F2**2._dl*(6._dl*F1*F2**2._dl*F3*(3._dl*F2**2._dl + F1*F3)*eft_cache%adotoa* &
        (-2._dl*F3**2._dl*F1dot + 3._dl*F2**2._dl*F3dot + 6._dl*F2*F3*(-F2dot + rho_d)) + &
        a*(6._dl*F2**4._dl*F3**3._dl*F1dot**2._dl + 2._dl*F1**3._dl*F3**2._dl*(-2._dl*F3*F2dot + F2*F3dot)**2._dl + &
        F1**2._dl*F2*F3*(-8._dl*F3**3._dl*F1dot*F2dot + F2**3._dl*(-9._dl*F3dot**2._dl + 6._dl*F3*F3dotdot) + &
        4._dl*F2*F3**2._dl*(F1dot*F3dot - F3*F1dotdot + 3._dl*(F2dot - rho_d)*rho_d) + &
        12._dl*F2**2._dl*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))) + &
        F1*(4._dl*F2**2._dl*F3**4._dl*F1dot**2._dl + 9._dl*F2**6._dl*(-3._dl*F3dot**2._dl + 2._dl*F3*F3dotdot) - &
        12._dl*F2**4._dl*F3**2._dl*(F3*F1dotdot + 3._dl*rho_d*(-F2dot + rho_d)) + &
        36._dl*F2**5._dl*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))))))/ &
        (a*(3._dl*F2**2._dl + F1*F3)*(3._dl*F2**2._dl + 2._dl*F1*F3)))/(16._dl*F1**2._dl*F2**4._dl*F3**2._dl)

    eft_cache%EFT_mu2 = ((-4._dl*F1*F2**2._dl*F3*(-2._dl*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2._dl)/ &
        (3._dl*F2**2._dl + F1*F3) + &
        (4._dl*F1*F3*(-2._dl*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2._dl)/(3._dl + (F1*F3)/F2**2._dl) - &
        (2._dl*F2**2._dl*F3*(6._dl*F1*F2**2._dl*F3*(3._dl*F2**2._dl + F1*F3)*eft_cache%adotoa* &
        (9._dl*F2**4._dl*F1dot + 6._dl*F1*F2**2._dl*F3*F1dot + &
        F1**2._dl*(-3._dl*F2**2._dl*F3dot + 2._dl*F3*(F3*F1dot + 3._dl*F2*(F2dot - rho_d)))) + &
        a*(-27._dl*F2**8._dl*F3*F1dot**2._dl - 2._dl*F1**5._dl*F3**2._dl*(-2._dl*F3*F2dot + F2*F3dot)**2._dl + &
        18._dl*F1*F2**6._dl*F3*(-2._dl*F3*F1dot**2._dl + 3._dl*F2**2._dl*F1dotdot) + &
        18._dl*F1**2._dl*F2**4._dl*F3* &
        (-(F1dot*(F3**2._dl*F1dot - 2._dl*F2*F3*F2dot + F2**2._dl*F3dot)) + 3._dl*F2**2._dl*F3*F1dotdot) + &
        2._dl*F1**3._dl*F2**2._dl*(-2._dl*F3**4._dl*F1dot**2._dl + 12._dl*F2*F3**3._dl*F1dot*F2dot + &
        9._dl*F2**4._dl*(F3dot**2._dl - F3*F3dotdot) + &
        6._dl*F2**2._dl*F3**2._dl*(-3._dl*F2dot**2._dl - F1dot*F3dot + 2._dl*F3*F1dotdot - 3._dl*F2dot*rho_d + &
        3._dl*rho_d**2._dl) + 18._dl*F2**3._dl*F3*(F3dot*rho_d + F3*(F2dotdot - rho_d_dot))) + &
        F1**4._dl*F2*F3*(8._dl*F3**3._dl*F1dot*F2dot + 3._dl*F2**3._dl*(F3dot**2._dl - 2._dl*F3*F3dotdot) + &
        4._dl*F2*F3**2._dl*(-(F1dot*F3dot) + F3*F1dotdot + 3._dl*(-2._dl*F2dot + rho_d)*(F2dot + rho_d)) + &
        12._dl*F2**2._dl*F3*(F3dot*(F2dot + rho_d) + F3*(F2dotdot - rho_d_dot))))))/ &
        (a*(3._dl*F2**2._dl + F1*F3)**2._dl*(3._dl*F2**2._dl + 2._dl*F1*F3)) + &
        (9._dl*F2**6._dl*(1._dl + 1._dl/(3._dl*sqrt(F2**4._dl/(3._dl*F2**2._dl + 2._dl*F1*F3)**2._dl)))* &
        (6._dl*F1*F2*F3*(3._dl*F2**2._dl + F1*F3)*eft_cache%adotoa* &
        (F2*F3*F1dot + F1*(F2*F3dot + 2._dl*F3*(-F2dot + rho_d))) + &
        a*(-3._dl*F2**4._dl*F3**2._dl*F1dot**2._dl + &
        2._dl*F1*F2**2._dl*F3**2._dl*(-(F3*F1dot**2._dl) + 3._dl*F2**2._dl*F1dotdot) + &
        2._dl*F1**3._dl*F3*(-2._dl*F2**2._dl*F3dot**2._dl + F2*F3*(F2*F3dotdot - 2._dl*F3dot*(-2._dl*F2dot + rho_d)) - &
        2._dl*F3**2*(F2dot**2._dl - F2dot*rho_d + rho_d**2._dl + F2*(F2dotdot - rho_d_dot))) + &
        F1**2*F2*(4._dl*F3**3._dl*F1dot*F2dot + F2**3._dl*(-9._dl*F3dot**2._dl + 6._dl*F3*F3dotdot) + &
        2._dl*F2*F3**2._dl*(-(F1dot*F3dot) + F3*F1dotdot + 6._dl*(F2dot - rho_d)*rho_d) + &
        12._dl*F2**2._dl*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))))))/ &
        (a*(3._dl*F2**2._dl + F1*F3)**2._dl*(3._dl*F2**2._dl + 2._dl*F1*F3)))/(16._dl*F1**2._dl*F2**4._dl*F3**2._dl)

    end subroutine EFTCAMBModelComputeMassFactors

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the quasi-static expression of mu (closed form valid up to Horndeski).
    !! The numerator and the denominator vanish together when pi does not couple to the metric
    !! (Omega' = c = Gamma2 = Gamma5 = 0): in that case the GR value mu = 1 is returned instead of 0/0.
    subroutine EFTCAMBModelComputeQSMu(self, a, k, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    real(dl), intent(in)                          :: k             !< the input scale.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: num, den, den_scale

    call EFTCAMBQSMuTerms(a, k, eft_par_cache, eft_cache, num, den, den_scale)
    if (EFTCAMBQSDegenerate(den, den_scale)) then
        eft_cache%QS_mu = 1._dl
    else
        eft_cache%QS_mu = num/den
    end if

    end subroutine EFTCAMBModelComputeQSMu

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the quasi-static expression of Sigma (lensing response).
    !! Sigma = mu (1 + gamma)/2 with mu the closed form quasi-static mu and gamma = Phi/Psi the
    !! quasi-static gravitational slip, computed from the same quasi-static coefficients used by
    !! EFTCAMBModelComputeQSMuSigma. A degenerate slip (0/0, no coupling of pi to the metric) is set to 1.
    subroutine EFTCAMBModelComputeQSSigma(self, a, k, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    real(dl), intent(in)                          :: k             !< the input scale.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: num, den, den_scale, mu, gamma
    real(dl) :: mu_num, mu_den, mu_den_scale, slip_num, slip_den, slip_den_scale

    ! quasi-static mu (closed form):
    call EFTCAMBQSMuTerms(a, k, eft_par_cache, eft_cache, num, den, den_scale)
    if (EFTCAMBQSDegenerate(den, den_scale)) then
        mu = 1._dl
    else
        mu = num/den
    end if
    ! quasi-static slip:
    call EFTCAMBQSGeneralTerms(a, k, eft_par_cache, eft_cache, mu_num, mu_den, mu_den_scale, &
        slip_num, slip_den, slip_den_scale)
    if (EFTCAMBQSDegenerate(slip_den, slip_den_scale)) then
        gamma = 1._dl
    else
        gamma = slip_num/slip_den
    end if

    eft_cache%QS_sigma = 0.5_dl*mu*(1._dl + gamma)

    end subroutine EFTCAMBModelComputeQSSigma

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the general quasi-static expressions of mu and Sigma.
    !! Degenerate (0/0) ratios are replaced by their GR value (mu = 1, gamma = 1).
    subroutine EFTCAMBModelComputeQSMuSigma(self, a, k, eft_par_cache, eft_cache) ! minsup

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    real(dl), intent(in)                          :: k             !< the input scale.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: mu_num, mu_den, mu_den_scale, slip_num, slip_den, slip_den_scale
    real(dl) :: gamma

    call EFTCAMBQSGeneralTerms(a, k, eft_par_cache, eft_cache, mu_num, mu_den, mu_den_scale, &
        slip_num, slip_den, slip_den_scale)

    if (EFTCAMBQSDegenerate(mu_den, mu_den_scale)) then
        eft_cache%QS_mu_gen = 1._dl
    else
        eft_cache%QS_mu_gen = mu_num/mu_den
    end if

    if (EFTCAMBQSDegenerate(slip_den, slip_den_scale)) then
        gamma = 1._dl
    else
        gamma = slip_num/slip_den
    end if

    eft_cache%QS_sigma_gen = 0.5_dl*eft_cache%QS_mu_gen*(1._dl + gamma)

    end subroutine EFTCAMBModelComputeQSMuSigma

    ! ---------------------------------------------------------------------------------------------
    !> Function that decides whether a quasi-static ratio is degenerate: its denominator vanishes
    !! within roundoff of the sum of the magnitudes of its terms (0/0 when all terms vanish).
    logical function EFTCAMBQSDegenerate(den, den_scale)

    implicit none

    real(dl), intent(in) :: den        !< the denominator of the ratio
    real(dl), intent(in) :: den_scale  !< the sum of the magnitudes of the terms of the denominator

    real(dl), parameter :: QS_relative_tolerance = 1.e-13_dl

    EFTCAMBQSDegenerate = .not. (abs(den) > QS_relative_tolerance*den_scale)

    end function EFTCAMBQSDegenerate

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes numerator and denominator of the closed form quasi-static mu,
    !! together with the sum of the magnitudes of the terms of the denominator.
    subroutine EFTCAMBQSMuTerms(a, k, eft_par_cache, eft_cache, num, den, den_scale)

    implicit none

    real(dl), intent(in)                          :: a             !< the input scale factor.
    real(dl), intent(in)                          :: k             !< the input scale.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    real(dl), intent(out) :: num        !< numerator of mu
    real(dl), intent(out) :: den        !< denominator of mu
    real(dl), intent(out) :: den_scale  !< sum of the magnitudes of the terms of the denominator

    real(dl) :: k2, H, Hd, Hdd, h0, opO, OP, c, G2, G2P, G5, G5P, X, Y, d1, d2, d3, d4, d5

    k2 = k**2
    H = eft_cache%adotoa
    Hd = eft_cache%Hdot
    Hdd = eft_cache%Hdotdot
    h0 = eft_par_cache%h0_Mpc
    opO = 1._dl + eft_cache%EFTOmegaV
    OP = eft_cache%EFTOmegaP
    c = eft_cache%EFTc
    G2 = eft_cache%EFTGamma2V
    G2P = eft_cache%EFTGamma2P
    G5 = eft_cache%EFTGamma5V
    G5P = eft_cache%EFTGamma5P

    num = 2._dl*(2._dl*opO*c*(k2 + 3._dl*H**2 - 3._dl*Hd) + &
        a*(-36._dl*H**4*opO*OP - 3._dl*H**3*h0*opO*(G2 + a*G2P) + &
        2._dl*k2*H**2*(2._dl*G5*OP + a*OP**2 + 2._dl*(opO + a*OP)*G5P) + &
        opO*(4._dl*k2*G5*Hd + 3._dl*h0*G2*Hdd) - &
        H*opO*(h0*(G2*(k2 + 3._dl*Hd) + a*(k2 - 3._dl*Hd)*G2P) - 12._dl*OP*Hdd)))

    X = opO + 2._dl*G5
    Y = a*h0*G2 + H*(-4._dl*G5 + a*OP)
    d1 = 2._dl*a*k2*H*X*OP*Y
    d2 = -k2*opO*Y**2
    d3 = 2._dl*k2*H*X*Y*(2._dl*G5 + a*(OP + 2._dl*G5P))
    d4 = 2._dl*k2*X**2*(2._dl*c - a*(-4._dl*G5*Hd + H*h0*(G2 + a*G2P) - 4._dl*H**2*G5P))
    d5 = 6._dl*X**2*(2._dl*c*(H**2 - Hd) - a*H*h0*(H**2 - Hd)*(3._dl*G2 + a*G2P) + &
        4._dl*a*H*OP*(-3._dl*H**3 + Hdd) + a*h0*G2*(2._dl*H**3 - 4._dl*H*Hd + Hdd))

    den = d1 + d2 + d3 + d4 + d5
    den_scale = abs(d1) + abs(d2) + abs(d3) + abs(d4) + abs(d5)

    end subroutine EFTCAMBQSMuTerms

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes numerators and denominators of the general quasi-static mu and
    !! gravitational slip gamma, together with the sums of the magnitudes of the denominator terms.
    subroutine EFTCAMBQSGeneralTerms(a, k, eft_par_cache, eft_cache, mu_num, mu_den, mu_den_scale, &
        slip_num, slip_den, slip_den_scale)

    implicit none

    real(dl), intent(in)                          :: a             !< the input scale factor.
    real(dl), intent(in)                          :: k             !< the input scale.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    real(dl), intent(out) :: mu_num, mu_den, mu_den_scale        !< general quasi-static mu
    real(dl), intent(out) :: slip_num, slip_den, slip_den_scale  !< quasi-static slip gamma

    real(dl) :: A1, A2, A3, B1, B2, B3, C1, C2, C3, M2, R0dot, r

    R0dot = (6._dl/a**3)*(eft_cache%Hdotdot - 2._dl*eft_cache%adotoa*(eft_cache%adotoa**2))

    A1 = 2._dl*(1._dl + eft_cache%EFTOmegaV) + 4._dl*eft_cache%EFTGamma5V
    A2 = -eft_cache%EFTOmegaP*eft_cache%adotoa - eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V + &
        2._dl*eft_cache%adotoa*eft_cache%EFTGamma4V/a + 4._dl*eft_cache%adotoa*eft_cache%EFTGamma5V/a
    A3 = -8._dl*eft_cache%EFTGamma6V

    B1 = -1._dl - 2._dl*eft_cache%EFTGamma5V/(1._dl + eft_cache%EFTOmegaV)
    B2 = 1._dl
    B3 = (-eft_cache%EFTOmegaP*eft_cache%adotoa + eft_cache%adotoa*(eft_cache%EFTGamma4V/a + eft_cache%EFTGamma4P)) &
        /(1._dl + eft_cache%EFTOmegaV)

    C1 = eft_cache%EFTOmegaP*eft_cache%adotoa + 2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V/a + eft_cache%EFTGamma5P)
    C2 = -eft_cache%EFTOmegaP*eft_cache%adotoa/2._dl - eft_cache%EFTGamma2V*eft_par_cache%h0_Mpc/2._dl - &
        (1._dl/(2._dl*a))*eft_cache%adotoa*(3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V - &
        4._dl*eft_cache%EFTGamma5V)
    C3 = eft_cache%EFTc/a**2 - eft_par_cache%h0_Mpc*eft_cache%adotoa*(eft_cache%EFTGamma2V/a + &
        eft_cache%EFTGamma2P)/2._dl &
        + eft_cache%EFTGamma3V*(k**2/2._dl - 3._dl*(eft_cache%Hdot - eft_cache%adotoa**2))/(a**2) &
        + eft_cache%EFTGamma4V*(k**2/2._dl - (eft_cache%Hdot - eft_cache%adotoa**2))/(a**2) &
        + 2._dl*(eft_cache%Hdot*eft_cache%EFTGamma5V/a**2 + eft_cache%adotoa**2*eft_cache%EFTGamma5P/a)

    M2 = eft_cache%EFTOmegaP*eft_cache%adotoa*R0dot/4._dl - &
        3._dl*eft_cache%EFTc*(eft_cache%Hdot - eft_cache%adotoa**2)/a**4 + &
        1.5_dl*(eft_par_cache%h0_Mpc/a**2)* &
        ((eft_cache%Hdot - eft_cache%adotoa**2)*eft_cache%adotoa*eft_cache%EFTGamma2P + &
        (eft_cache%Hdotdot - eft_cache%adotoa*eft_cache%Hdot - eft_cache%adotoa**3)*eft_cache%EFTGamma2V/a) + &
        1.5_dl*((eft_cache%Hdot - eft_cache%adotoa**2)**2/a**4)*(eft_cache%EFTGamma4V + 3._dl*eft_cache%EFTGamma3V)

    r = a**2/k**2

    mu_num = 2._dl*(B2*C3 - C1*B3 + r*B2*M2)
    mu_den = A1*(B3*C2 - B1*C3 - B1*M2*r) + A2*(B1*C1 - B2*C2) + A3*(B2*C3 - B3*C1 + B2*M2*r)
    mu_den_scale = abs(A1)*(abs(B3*C2) + abs(B1*C3) + abs(B1*M2*r)) + abs(A2)*(abs(B1*C1) + abs(B2*C2)) + &
        abs(A3)*(abs(B2*C3) + abs(B3*C1) + abs(B2*M2*r))

    slip_num = B3*C2 - B1*C3 - r*B1*M2
    slip_den = B2*C3 - C1*B3 + r*B2*M2
    slip_den_scale = abs(B2*C3) + abs(C1*B3) + abs(r*B2*M2)

    end subroutine EFTCAMBQSGeneralTerms

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBModelAdditionalModelStability(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.
    logical :: EFTCAMBModelAdditionalModelStability

    EFTCAMBModelAdditionalModelStability = .true.

    end function EFTCAMBModelAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the LHS of positivity bounds. For details refer to the master thesis of Dani de Boe.
    subroutine EFTCAMBModelPositivityBounds(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: a2, one_plus_omega, adotoa2, adotoa02, NumFactor
    real(dl) :: D1, N1, N2, N3, N4, N5, N6
    real(dl) :: Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, tach_pos_num

    a2 = a**2
    one_plus_omega = 1._dl + eft_cache%EFTOmegaV
    adotoa2 = eft_cache%adotoa**2
    adotoa02 = eft_par_cache%h0_Mpc**2
    NumFactor = 2*a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V + &
        eft_cache%adotoa*(4*eft_cache%EFTGamma3V + a*eft_cache%EFTGamma3P)

    D1 = (-2*eft_cache%EFTc + 3*a*eft_cache%adotoa*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V + &
        (adotoa2 + eft_cache%Hdot*2.)*eft_cache%EFTGamma3V + 2*a*adotoa2*eft_cache%EFTGamma3P)

    N1 = D1 + a*(eft_cache%Hdot*eft_cache%EFTOmegaP + a*adotoa2*eft_cache%EFTOmegaPP)

    N2 = 4*a2*adotoa02*eft_cache%EFTGamma1V + &
        a*eft_par_cache%h0_Mpc*eft_cache%adotoa*(-3*eft_cache%EFTGamma2V + a*eft_cache%EFTGamma2P) + &
        eft_cache%Hdot*(4*eft_cache%EFTGamma3V - a*eft_cache%EFTGamma3P) + &
        adotoa2*(-10*eft_cache%EFTGamma3V + a*(eft_cache%EFTGamma3P - a*eft_cache%EFTGamma3PP))

    N3 = -2*eft_cache%EFTcdot + 3*a*eft_par_cache%h0_Mpc* &
        (eft_cache%EFTGamma2V*eft_cache%Hdot + adotoa2*(-eft_cache%EFTGamma2V + a*eft_cache%EFTGamma2P)) - &
        2*eft_cache%Hdot*eft_cache%adotoa*(eft_cache%EFTGamma3V - 3*a*eft_cache%EFTGamma3P) + &
        2*eft_cache%EFTGamma3V*eft_cache%Hdotdot - &
        eft_cache%adotoa**3.*(2*eft_cache%EFTGamma3V + a*eft_cache%EFTGamma3P - 2*a2*eft_cache%EFTGamma3PP)

    N4 = -2*eft_cache%EFTLambdadot*eft_cache%adotoa + 2*eft_cache%EFTLambdadotdot - &
        a2*eft_cache%adotoa*eft_par_cache%h0_Mpc*eft_cache%Hdot*(eft_cache%EFTGamma2P - 3*a*eft_cache%EFTGamma2PP) - &
        eft_cache%adotoa*eft_cache%Hdotdot* &
        (-2*eft_cache%EFTGamma3V + 5*a*eft_cache%EFTGamma3P + 4*a2*eft_cache%EFTGamma3PP) - &
        2*eft_cache%EFTGamma3V*(2*eft_cache%Hdot**2. + eft_cache%Hdotdotdot)

    N5 = a2*eft_par_cache%h0_Mpc*eft_cache%Hdotdot*eft_cache%EFTGamma2P - &
        a*eft_cache%Hdot**2.*(5*eft_cache%EFTGamma3P + 3*a*eft_cache%EFTGamma3PP) - &
        a*eft_cache%Hdotdotdot*eft_cache%EFTGamma3P + &
        a2**2.*eft_par_cache%h0_Mpc*eft_cache%adotoa**3.*eft_cache%EFTGamma2PPP + &
        36*adotoa2*eft_cache%Hdot*eft_cache%EFTGamma3V - &
        2*adotoa2*a2*eft_cache%Hdot*(4*eft_cache%EFTGamma3PP + 3*a*eft_cache%EFTGamma3PPP)

    N6 = adotoa2**2.*(24*eft_cache%EFTGamma3V - 12*a*eft_cache%EFTGamma3P + 2*a**3.*eft_cache%EFTGamma3PPP + &
        a2**2.*eft_cache%EFTGamma3PPPP)

    ! ghost-positivity condition:
    eft_cache%ghostpos = D1

    ! tachyonic-positivity condition:
    Q1 = (-eft_cache%adotoa*eft_cache%EFTLambdadot + eft_cache%EFTLambdadotdot)*a2
    Q2 = eft_par_cache%h0_Mpc*(eft_cache%Hdotdot - eft_cache%adotoa*eft_cache%Hdot)*a2*eft_cache%EFTGamma2P/2._dl
    Q3 = 3._dl*eft_par_cache%h0_Mpc*eft_cache%adotoa*eft_cache%Hdot*(a**3)*eft_cache%EFTGamma2PP/2._dl
    Q4 = eft_par_cache%h0_Mpc*(eft_cache%adotoa**3)*(a2**2)*eft_cache%EFTGamma2PPP/2._dl
    Q5 = (-eft_cache%Hdotdotdot - 12._dl*eft_cache%adotoa**4 + 18._dl*eft_cache%adotoa**2*eft_cache%Hdot - &
        2._dl*eft_cache%Hdot**2 + eft_cache%adotoa*eft_cache%Hdotdot)*eft_cache%EFTGamma3V
    Q6 = (a/2._dl)*(-eft_cache%Hdotdotdot + 12._dl*eft_cache%adotoa**4 - 5._dl*eft_cache%Hdot**2 - &
        5._dl*eft_cache%adotoa*eft_cache%Hdotdot)*eft_cache%EFTGamma3P
    Q7 = -(a2/2._dl)*(3._dl*eft_cache%Hdot**2 + 8._dl*eft_cache%adotoa**2*eft_cache%Hdot + &
        4._dl*eft_cache%Hdotdot*eft_cache%adotoa)*eft_cache%EFTGamma3PP
    Q8 = -(a**3)*eft_cache%adotoa**2*(3._dl*eft_cache%Hdot + eft_cache%adotoa**2)*eft_cache%EFTGamma3PPP
    Q9 = -a2**2*eft_cache%adotoa**4*eft_cache%EFTGamma3PPPP/(2._dl)
    tach_pos_num = Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9

    ! compute the LHS of positivity bounds:
    if (D1 /= 0) then

        eft_cache%posbound1 = -(4*a2*one_plus_omega*adotoa02*eft_cache%EFTGamma2V**2. + &
            4*eft_cache%EFTGamma3V**3.*(adotoa2 + 2*eft_cache%Hdot) + &
            8*a*adotoa2*one_plus_omega*eft_cache%EFTGamma3V*eft_cache%EFTGamma3P + &
            a2*adotoa2*one_plus_omega*eft_cache%EFTGamma3P**2. + &
            4*a*eft_cache%adotoa*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V* &
            (4*one_plus_omega*eft_cache%EFTGamma3V + 3*eft_cache%EFTGamma3V**2. + &
            a*one_plus_omega*eft_cache%EFTGamma3P) + &
            8*eft_cache%EFTGamma3V**2.*(-eft_cache%EFTc + adotoa2*(2*one_plus_omega + a*eft_cache%EFTGamma3P)))/D1/8* &
            one_plus_omega

        eft_cache%posbound2 = 0.5*(2*eft_cache%EFTGamma3V*N1/one_plus_omega + &
            2*a*(eft_cache%Hdot*eft_cache%EFTGamma3P + a*adotoa2*eft_cache%EFTGamma3PP) + N2 + NumFactor*-N3/D1 + &
            NumFactor**2.*(-N4 - N5 + N6)/D1**2./8)
        eft_cache%tachpos = -tach_pos_num/D1
    else
        eft_cache%posbound1 = 0
        eft_cache%posbound2 = 0
        eft_cache%tachpos = 0
    end if
    end subroutine EFTCAMBModelPositivityBounds

    function EFTCAMBModelComputeDeltaTau(self, a1, a2)

    implicit none

    class(EFTCAMB_model) :: self          !< the base class
    real(dl), intent(in) :: a1            !< the input scale factor.
    real(dl), intent(in) :: a2            !< the input scale factor.

    real(dl)             :: EFTCAMBModelComputeDeltaTau !< the output tau

    write(*, *) "EFTCAMB_use_background = True"
    write(*, *) "EFTCAMBModelComputeDeltaTau currently only supported by the Horndeski module (EFTflag=5)"
    call MpiStop("EFTCAMB error")

    end function EFTCAMBModelComputeDeltaTau

    function EFTCAMBModelComputeDeltaCosmicT(self, a1, a2)

    implicit none

    class(EFTCAMB_model) :: self          !< the base class
    real(dl), intent(in) :: a1            !< the input scale factor.
    real(dl), intent(in) :: a2            !< the input scale factor.

    real(dl)             :: EFTCAMBModelComputeDeltaCosmicT !< the output t

    write(*, *) "EFTCAMB_use_background = True"
    write(*, *) "EFTCAMBModelComputeDeltaCosmicT currently only supported by the Horndeski module (EFTflag=5)"
    call MpiStop("EFTCAMB error")

    end function EFTCAMBModelComputeDeltaCosmicT

    function EFTCAMBModelComputeSoundHorizon(self, a)

    implicit none

    class(EFTCAMB_model) :: self          !< the base class
    real(dl), intent(in) :: a             !< the input scale factor.

    real(dl)             :: EFTCAMBModelComputeSoundHorizon !< the output r_s

    write(*, *) "EFTCAMB_use_background = True"
    write(*, *) "EFTCAMBModelComputeSoundHorizon currently only supported by the Horndeski module (EFTflag=5)"
    call MpiStop("EFTCAMB error")

    end function EFTCAMBModelComputeSoundHorizon

    end module EFTCAMB_abstract_model

! ----------------------------------------------------------------------------------------

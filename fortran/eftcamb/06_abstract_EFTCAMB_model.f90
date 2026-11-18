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
    use config, only: GlobalError, error_unsupported_params
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
        !< effective CPL (w0, wa) of the model background, computed by EFTCAMB_compute_effective_w0wa
        !! (09_EFTCAMB_stability.f90) after the background initialization; used when the model does not
        !! provide its own (w0, wa).
        logical                       :: background_w0wa_computed = .false.
        real(dl)                      :: w0_background = -1._dl
        real(dl)                      :: wa_background = 0._dl

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
    !< subroutine that fills the alpha-basis outputs (alphaB/K/M/T, their derivatives, Meff2) from the EFT functions.
    procedure :: compute_alpha_functions => EFTCAMBModelComputeAlphaFunctions
    !< subroutine that computes the Einstein equations factors. For details refer to the numerical notes.
    procedure :: compute_Einstein_factors => EFTCAMBModelComputeEinsteinFactors
    !< subroutine that computes the pi field equations factors. For details refer to the numerical notes.
    procedure :: compute_pi_factors => EFTCAMBModelComputePiFactors
    !< subroutine that computes the factors for the tensor propagation equation. For details refer to the numerical notes.
    procedure :: compute_tensor_factors => EFTCAMBModelComputeTensorFactors
    !< subroutine that computes the kinetic and gradient terms. For details refer to the numerical notes.
    procedure :: compute_stability_factors => EFTCAMBModelComputeStabilityFactors
    !< subroutine that computes the k -> 0 mass eigenvalues of the scalar sector (arXiv:1609.03599 App. C, corrected).
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
    !< subroutine that tests the model's own scalar ghost and gradient conditions, for models beyond GLPV.
    procedure :: beyond_GLPV_stability => EFTCAMBModelBeyondGLPVStability
    !< function that tests the phenomenological prior w_DE <= -1/3 (EFT_wDE_prior).
    procedure :: wDE_prior => EFTCAMBModelwDEPrior
    !< function that returns the scale factor where the model background solution starts (0: everywhere).
    procedure :: background_start => EFTCAMBModelBackgroundStart
    !< function that reads model_background_a_final (the end of a model background table) and returns its log.
    procedure :: read_background_x_final => EFTCAMBModelReadBackgroundXFinal
    !< subroutine that checks that the pi evolution is regular after the perturbation turn-on (message: '' if it is).
    procedure :: pi_evolution_check => EFTCAMBModelPiEvolutionCheck
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
    class(TIniFile)      :: Ini       !< Input ini file
    integer              :: eft_error !< error code: 0 all fine, 1 initialization failed
    end subroutine EFTCAMBModelReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBModelAllocateModelSelection(self, Ini, eft_error)
    use IniObjects
    import EFTCAMB_model
    implicit none
    class(EFTCAMB_model) :: self      !< the base class
    class(TIniFile)      :: Ini       !< Input ini file for additional initialization
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
    class(TIniFile)      :: Ini       !< Input ini file
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
    !> Subroutine that fills the alpha-basis outputs from the EFT functions (outputs only: get_eft_functions and the
    !! evolution). Bellini-Sawicki arXiv:1404.3713 mapped as in Gleyzes et al. arXiv:1411.3712, in the cache
    !! convention of the models that compute them themselves (RPH and Horndeski override this subroutine; these
    !! formulas reproduce their values to roundoff and 1e-8):
    !!   Meff2 = M*^2/m0^2 = 1 + Omega - gamma4,  alphaT = gamma4/Meff2,  alphaM = a (Omega' - gamma4')/Meff2,
    !!   alphaB = (a Hc Omega' + a H0 gamma2)/(2 Hc Meff2)  (= -alpha_B^BS/2, the EFTCAMB/RPH convention),
    !!   alphaK = (2 c + 4 a^2 H0^2 gamma1)/(Hc^2 Meff2)  (cache c = a^2 c),
    !! and their conformal time derivatives a Hc d/da. gamma4'' is the model's, or -gamma3'' on the Horndeski line
    !! when the model leaves it 0.
    subroutine EFTCAMBModelComputeAlphaFunctions(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout)  :: eft_cache

    real(dl) :: Hc, H0, g3, g4, g4pp, M, Mp, Mpp, N, Np, dHda, K, Kdot

    Hc = eft_cache%adotoa
    H0 = eft_par_cache%h0_Mpc
    g3 = eft_cache%EFTGamma3V
    g4 = eft_cache%EFTGamma4V
    g4pp = eft_cache%EFTGamma4PP
    if (g4pp == 0._dl .and. abs(g3 + g4) <= 1.e-12_dl*(abs(g3) + abs(g4))) g4pp = -eft_cache%EFTGamma3PP
    M = 1._dl + eft_cache%EFTOmegaV - g4
    if (M == 0._dl .or. Hc == 0._dl) return
    Mp = eft_cache%EFTOmegaP - eft_cache%EFTGamma4P
    Mpp = eft_cache%EFTOmegaPP - g4pp
    eft_cache%Meff2 = M
    eft_cache%alphaT = g4/M
    eft_cache%alphaTdot = a*Hc*(eft_cache%EFTGamma4P/M - g4*Mp/M**2)
    eft_cache%alphaM = a*Mp/M
    eft_cache%alphaMdot = a*Hc*((Mp + a*Mpp)/M - a*Mp**2/M**2)
    N = a*Hc*eft_cache%EFTOmegaP + a*H0*eft_cache%EFTGamma2V
    dHda = eft_cache%Hdot/(a*Hc)
    Np = Hc*eft_cache%EFTOmegaP + a*dHda*eft_cache%EFTOmegaP + a*Hc*eft_cache%EFTOmegaPP + H0*eft_cache%EFTGamma2V + &
        a*H0*eft_cache%EFTGamma2P
    eft_cache%alphaB = N/(2._dl*Hc*M)
    eft_cache%alphaBdot = a*Hc*(Np/(2._dl*Hc*M) - N*dHda/(2._dl*Hc**2*M) - N*Mp/(2._dl*Hc*M**2))
    K = 2._dl*eft_cache%EFTc + 4._dl*a**2*H0**2*eft_cache%EFTGamma1V
    Kdot = 4._dl*Hc*eft_cache%EFTc + 2._dl*eft_cache%EFTcdot + 8._dl*a**2*Hc*H0**2*eft_cache%EFTGamma1V + &
        4._dl*a**3*Hc*H0**2*eft_cache%EFTGamma1P
    eft_cache%alphaK = K/(Hc**2*M)
    eft_cache%alphaKdot = Kdot/(Hc**2*M) - 2._dl*K*eft_cache%Hdot/(Hc**3*M) - K*a*Hc*Mp/(Hc**2*M**2)

    end subroutine EFTCAMBModelComputeAlphaFunctions

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
    ! The terms above are the bare form: every matter species in the alpha_H^2 mixing, as for fluids with pressure.
    ! Pressureless dust dresses it (decision D3): in GLPV the DE mode propagates with c_s^2 + lambda_d,
    ! lambda_d = (2 G5 - F1)^2 rho_dust/(F2^2 L_zz) (Kase & Tsujikawa arXiv:1802.02728; measured on the evolved pi:
    ! c^2 = 0.381 vs 0.384), so the gradient condition is G + (2 G5 - F1)^2 rho_dust/F2^2 >= 0 with
    ! G5 = 1 + Omega + 2 gamma5, F1 = 2 (1 + Omega) + 3 gamma3 + gamma4. The dressing vanishes in Horndeski (2 G5 = F1).
    eft_cache%EFT_gradient = eft_cache%EFT_gradient + 4.5_dl*(4._dl*eft_cache%EFTGamma5V - &
        3._dl*eft_cache%EFTGamma3V - eft_cache%EFTGamma4V)**2*(eft_cache%grhob_t + eft_cache%grhoc_t)

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
        + eft_cache%EFTOmegaV)* (4._dl + a*eft_cache%EFTOmegaP + 4._dl*eft_cache%EFTOmegaV)*eft_cache%Hdot))) + &
        4.5_dl*(4._dl*eft_cache%EFTGamma5V - 3._dl*eft_cache%EFTGamma3V - eft_cache%EFTGamma4V)**2* &
        (eft_cache%grhob_t + eft_cache%grhoc_t)

    end subroutine EFTCAMBStabilityScales

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the k -> 0 mass eigenvalues of the scalar sector, EFT_mu1 and EFT_mu2 (cosmic time,
    !! Mpc^-2), from the closed form of arXiv:1609.03599 App. C (pi and dust at k -> 0, after diagonalising the kinetic
    !! matrix). In exact arithmetic the closed form is right: it gives V''/2 for canonical quintessence and half the
    !! f(R) scalaron mass, with its sign. Three corrections make it usable (the original body is kept below, commented):
    !!  - EFTcPP = d^2 EFTc/da^2 (the original local expression was too large by about 1/adotoa);
    !!  - the closed form assumes a dust-only background. It gets the rho + p of all matter as rho_d, with dust
    !!    continuity, and the Hdotdot of the background in which that dust obeys the c equation and its time
    !!    derivative: Hdotdot + a^2 dP/dtau/(2 (1 + Omega) + a dOmega/da). Otherwise the non-dust residual is divided by
    !!    Ktilde = 3 F2^2 + F1 F3, which vanishes near GR (radiation at 1e-4 of the matter density gives 200 V'' for
    !!    quintessence). The correction leaves an error of order rho_rad/rho_m (1e-4 at 1e-4, 13% at 1e-1);
    !!  - Ktilde cancels near GR, so the closed form is evaluated in quadruple precision (double where the compiler
    !!    has none) and only where Ktilde/(3 F2^2) > 100 sqrt(epsilon), which keeps the roundoff error below
    !!    1e-3 H^2. Elsewhere (and at a NaN) resolved = .false. and EFT_mu1 = EFT_mu2 = 0.
    !! The outputs are twice the eigenvalues of the paper, so that EFT_mu1 is the mass squared (V'' for canonical
    !! quintessence, the f(R) scalaron mass ((1 + f_R)/f_RR - R)/3): a mode grows as exp(int sqrt(-mu) dt) where
    !! mu < 0. EFT_mu2 is the matter (dust) eigenvalue.
    subroutine EFTCAMBModelComputeMassFactors(self, a, eft_par_cache, eft_cache, resolved)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    logical, intent(out), optional                :: resolved      !< whether EFT_mu1, EFT_mu2 are resolved.

    ! quadruple precision where the compiler has it:
    integer, parameter :: qp = merge(selected_real_kind(30), dl, selected_real_kind(30) > 0)
    real(qp) :: F1, F2, F3, F1dot, F2dot, rho_d, rho_d_dot, F1dotdot, F2dotdot, F3dot, F3dotdot
    real(qp) :: EFTcP, EFTcPP, Ktilde, denom, mu1, mu2
    real(qp) :: h0_Mpc, Hdotdot, EFTGamma1P, EFTGamma1PP, EFTGamma1V, EFTGamma2P, EFTGamma2PP, EFTGamma2V, &
        EFTGamma3P, EFTGamma3PP, EFTGamma3V, EFTGamma4P, EFTGamma4PP, EFTGamma4V, EFTOmegaP, EFTOmegaPP, &
        EFTOmegaPPP, EFTOmegaV, EFTc, EFTcdot, EFTcdotdot, Hdot, adotoa

    eft_cache%EFT_mu1 = 0._dl
    eft_cache%EFT_mu2 = 0._dl
    if (present(resolved)) resolved = .false.

    ! inputs:
    h0_Mpc = eft_par_cache%h0_Mpc
    EFTGamma1P = eft_cache%EFTGamma1P
    EFTGamma1PP = eft_cache%EFTGamma1PP
    EFTGamma1V = eft_cache%EFTGamma1V
    EFTGamma2P = eft_cache%EFTGamma2P
    EFTGamma2PP = eft_cache%EFTGamma2PP
    EFTGamma2V = eft_cache%EFTGamma2V
    EFTGamma3P = eft_cache%EFTGamma3P
    EFTGamma3PP = eft_cache%EFTGamma3PP
    EFTGamma3V = eft_cache%EFTGamma3V
    EFTGamma4P = eft_cache%EFTGamma4P
    EFTGamma4PP = eft_cache%EFTGamma4PP
    EFTGamma4V = eft_cache%EFTGamma4V
    EFTOmegaP = eft_cache%EFTOmegaP
    EFTOmegaPP = eft_cache%EFTOmegaPP
    EFTOmegaPPP = eft_cache%EFTOmegaPPP
    EFTOmegaV = eft_cache%EFTOmegaV
    EFTc = eft_cache%EFTc
    EFTcdot = eft_cache%EFTcdot
    EFTcdotdot = eft_cache%EFTcdotdot
    Hdot = eft_cache%Hdot
    adotoa = eft_cache%adotoa

    ! dust-equivalent background: the matter rho + p as dust, and Hdotdot shifted so that the c equation and its time
    ! derivative hold with dust continuity:
    denom = 2._qp*(1._qp + EFTOmegaV) + a*EFTOmegaP
    if (denom == 0._qp) return
    rho_d = (real(eft_cache%grhom_t, qp) + real(eft_cache%gpresm_t, qp))/a**2
    rho_d_dot = -3._qp*adotoa/a*rho_d
    Hdotdot = real(eft_cache%Hdotdot, qp) + real(eft_cache%gpresdotm_t, qp)/denom

    ! d EFTc/da and d^2 EFTc/da^2:
    EFTcP = EFTcdot/(a*adotoa) + 2._qp*EFTc/a
    EFTcPP = (EFTcdotdot/adotoa + (3._qp - Hdot/adotoa**2)*EFTcdot + 2._qp*adotoa*EFTc)/(a**2*adotoa)

    ! the closed form of arXiv:1609.03599 App. C:
    F1 = 2._dl*(1._dl + EFTOmegaV) + 3._dl*EFTGamma3V + EFTGamma4V

    F2 = h0_Mpc*EFTGamma2V + &
        (adotoa*(2._dl*(1._dl + EFTOmegaV) + 3._dl*EFTGamma3V + EFTGamma4V))/ &
        a + adotoa*EFTOmegaP

    F3 = (2._dl*EFTc)/a**2 + 4._dl*h0_Mpc**2*EFTGamma1V - &
        (6._dl*h0_Mpc*adotoa*EFTGamma2V)/a - &
        (3._dl*adotoa**2* &
        (2._dl*(1._dl + EFTOmegaV) + 3._dl*EFTGamma3V + EFTGamma4V))/a**2 - &
        (6._dl*adotoa**2*EFTOmegaP)/a

    F1dot = adotoa*(2._dl*EFTOmegaP + 3._dl*EFTGamma3P + EFTGamma4P)

    F2dot = ((2._dl + 2._dl*EFTOmegaV + 3._dl*EFTGamma3V + EFTGamma4V + &
        a*EFTOmegaP)*Hdot + &
        a**2*h0_Mpc**2*adotoa*EFTGamma2P + &
        adotoa**2*(-2._dl - 2._dl*EFTOmegaV - 3._dl*EFTGamma3V - &
        EFTGamma4V + 2._dl*a*EFTOmegaP + 3._dl*a*EFTGamma3P + a*EFTGamma4P + &
        a**2*EFTOmegaPP))/a**2

    F3dot = -((4._dl*EFTc*adotoa + &
        6._dl*a*h0_Mpc*EFTGamma2V*Hdot - &
        2._dl*adotoa*(a*EFTcP - &
        3._dl*(2._dl + 2._dl*EFTOmegaV + 3._dl*EFTGamma3V + EFTGamma4V + &
        2._dl*a*EFTOmegaP)*Hdot + 2._dl*a**3*h0_Mpc**2*EFTGamma1P) + &
        6._dl*a*h0_Mpc*adotoa**2*(-EFTGamma2V + a*EFTGamma2P) - &
        3._dl*adotoa**3* &
        (4._dl + 4._dl*EFTOmegaV + 6._dl*EFTGamma3V + 2._dl*EFTGamma4V - &
        3._dl*a*EFTGamma3P - a*EFTGamma4P - 2._dl*a**2*EFTOmegaPP))/a**3)

    F1dotdot = (2._dl*EFTOmegaP*Hdot + &
        Hdot*(3._dl*EFTGamma3P + EFTGamma4P) + &
        a*adotoa**2*(2._dl*EFTOmegaPP + 3._dl*EFTGamma3PP + EFTGamma4PP))/a

    F2dotdot = (a**2*h0_Mpc**2*Hdot*EFTGamma2P + &
        adotoa*Hdot* &
        (-8._dl - 8._dl*EFTOmegaV - 12._dl*EFTGamma3V - 4._dl*EFTGamma4V + &
        5._dl*a*EFTOmegaP + 9._dl*a*EFTGamma3P + 3._dl*a*EFTGamma4P + &
        3._dl*a**2*EFTOmegaPP) + &
        (2._dl + 2._dl*EFTOmegaV + 3._dl*EFTGamma3V + EFTGamma4V + &
        a*EFTOmegaP)*Hdotdot + &
        a**3*h0_Mpc**2*adotoa**2*EFTGamma2PP + &
        adotoa**3*(4._dl + 4._dl*EFTOmegaV + 6._dl*EFTGamma3V + &
        2._dl*EFTGamma4V - 4._dl*a*EFTOmegaP - 6._dl*a*EFTGamma3P - &
        2._dl*a*EFTGamma4P + 2._dl*a**2*EFTOmegaPP + 3._dl*a**2*EFTGamma3PP + &
        a**2*EFTGamma4PP + a**3*EFTOmegaPPP))/a**3

    F3dotdot = (4._dl*EFTc*(3._dl*adotoa**2 - Hdot) + &
        2._dl*(a*EFTcP*Hdot - &
        3._dl*(2._dl + 2._dl*EFTOmegaV + 3._dl*EFTGamma3V + EFTGamma4V + &
        2*a*EFTOmegaP)*Hdot**2 + &
        2._dl*a**3*h0_Mpc**2*Hdot*EFTGamma1P - &
        3._dl*a*h0_Mpc*EFTGamma2V*Hdotdot) - &
        6._dl*adotoa*(a*h0_Mpc*Hdot* &
        (-4._dl*EFTGamma2V + 3._dl*a*EFTGamma2P) + &
        (2._dl + 2._dl*EFTOmegaV + 3._dl*EFTGamma3V + EFTGamma4V + &
        2._dl*a*EFTOmegaP)*Hdotdot) + &
        adotoa**2*(-8._dl*a*EFTcP + &
        3._dl*Hdot*(24._dl + 24._dl*EFTOmegaV + 36._dl*EFTGamma3V + &
        12._dl*EFTGamma4V + 4._dl*a*EFTOmegaP - 15._dl*a*EFTGamma3P - &
        5._dl*a*EFTGamma4P - 10._dl*a**2*EFTOmegaPP) + &
        2._dl*a**2*(EFTcPP + 2._dl*a**2*h0_Mpc**2*EFTGamma1PP)) - &
        6._dl*a*h0_Mpc*adotoa**3* &
        (2*EFTGamma2V + a*(-2._dl*EFTGamma2P + a*EFTGamma2PP)) - &
        3._dl*adotoa**4* &
        (12._dl + 12._dl*EFTOmegaV + 18._dl*EFTGamma3V + 6._dl*EFTGamma4V - &
        4._dl*a*EFTOmegaP - 12._dl*a*EFTGamma3P - 4._dl*a*EFTGamma4P - &
        2._dl*a**2*EFTOmegaPP + 3._dl*a**2*EFTGamma3PP + a**2*EFTGamma4PP + &
        2._dl*a**3*EFTOmegaPPP))/a**4


    ! Ktilde = 3 F2^2 + F1 F3 divides the eigenvalues: resolved only well above roundoff.
    Ktilde = 3._qp*F2**2 + F1*F3
    if (.not. (abs(Ktilde) > 100._qp*sqrt(epsilon(1._qp))*3._qp*F2**2)) return

    mu1 = ((-4._dl*F1*F2**2*F3*(-2._dl*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2)/ &
        (3._dl*F2**2 + F1*F3) + &
        (4._dl*F1*F3*(-2._dl*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2)/(3._dl + (F1*F3)/F2**2) - &
        (9._dl*F2**6*(1._dl + 1._dl/(3._dl*sqrt(F2**4/(3._dl*F2**2 + 2._dl*F1*F3)**2)))* &
        (6._dl*F1*F2*F3*(3._dl*F2**2 + F1*F3)*adotoa* &
        (F2*F3*F1dot + F1*(F2*F3dot + 2._dl*F3*(-F2dot + rho_d))) + &
        a*(-3._dl*F2**4*F3**2*F1dot**2 + &
        2._dl*F1*F2**2*F3**2*(-(F3*F1dot**2) + 3._dl*F2**2*F1dotdot) + &
        2._dl*F1**3*F3*(-2._dl*F2**2*F3dot**2 + F2*F3*(F2*F3dotdot - 2._dl*F3dot*(-2._dl*F2dot + rho_d)) - &
        2._dl*F3**2*(F2dot**2 - F2dot*rho_d + rho_d**2 + F2*(F2dotdot - rho_d_dot))) + &
        F1**2*F2*(4._dl*F3**3*F1dot*F2dot + F2**3*(-9._dl*F3dot**2 + 6._dl*F3*F3dotdot) + &
        2._dl*F2*F3**2*(-(F1dot*F3dot) + F3*F1dotdot + 6._dl*(F2dot - rho_d)*rho_d) + &
        12._dl*F2**2*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))))))/ &
        (a*(3._dl*F2**2 + F1*F3)**2*(3._dl*F2**2 + 2._dl*F1*F3)) + &
        (2._dl*F1*F2**2*(6._dl*F1*F2**2*F3*(3._dl*F2**2 + F1*F3)*adotoa* &
        (-2._dl*F3**2*F1dot + 3._dl*F2**2*F3dot + 6._dl*F2*F3*(-F2dot + rho_d)) + &
        a*(6._dl*F2**4*F3**3*F1dot**2 + 2._dl*F1**3*F3**2*(-2._dl*F3*F2dot + F2*F3dot)**2 + &
        F1**2*F2*F3*(-8._dl*F3**3*F1dot*F2dot + F2**3*(-9._dl*F3dot**2 + 6._dl*F3*F3dotdot) + &
        4._dl*F2*F3**2*(F1dot*F3dot - F3*F1dotdot + 3._dl*(F2dot - rho_d)*rho_d) + &
        12._dl*F2**2*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))) + &
        F1*(4._dl*F2**2*F3**4*F1dot**2 + 9._dl*F2**6*(-3._dl*F3dot**2 + 2._dl*F3*F3dotdot) - &
        12._dl*F2**4*F3**2*(F3*F1dotdot + 3._dl*rho_d*(-F2dot + rho_d)) + &
        36._dl*F2**5*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))))))/ &
        (a*(3._dl*F2**2 + F1*F3)*(3._dl*F2**2 + 2._dl*F1*F3)))/(16._dl*F1**2*F2**4*F3**2)

    mu2 = ((-4._dl*F1*F2**2*F3*(-2._dl*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2)/ &
        (3._dl*F2**2 + F1*F3) + &
        (4._dl*F1*F3*(-2._dl*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2)/(3._dl + (F1*F3)/F2**2) - &
        (2._dl*F2**2*F3*(6._dl*F1*F2**2*F3*(3._dl*F2**2 + F1*F3)*adotoa* &
        (9._dl*F2**4*F1dot + 6._dl*F1*F2**2*F3*F1dot + &
        F1**2*(-3._dl*F2**2*F3dot + 2._dl*F3*(F3*F1dot + 3._dl*F2*(F2dot - rho_d)))) + &
        a*(-27._dl*F2**8*F3*F1dot**2 - 2._dl*F1**5*F3**2*(-2._dl*F3*F2dot + F2*F3dot)**2 + &
        18._dl*F1*F2**6*F3*(-2._dl*F3*F1dot**2 + 3._dl*F2**2*F1dotdot) + &
        18._dl*F1**2*F2**4*F3* &
        (-(F1dot*(F3**2*F1dot - 2._dl*F2*F3*F2dot + F2**2*F3dot)) + 3._dl*F2**2*F3*F1dotdot) + &
        2._dl*F1**3*F2**2*(-2._dl*F3**4*F1dot**2 + 12._dl*F2*F3**3*F1dot*F2dot + &
        9._dl*F2**4*(F3dot**2 - F3*F3dotdot) + &
        6._dl*F2**2*F3**2*(-3._dl*F2dot**2 - F1dot*F3dot + 2._dl*F3*F1dotdot - 3._dl*F2dot*rho_d + &
        3._dl*rho_d**2) + 18._dl*F2**3*F3*(F3dot*rho_d + F3*(F2dotdot - rho_d_dot))) + &
        F1**4*F2*F3*(8._dl*F3**3*F1dot*F2dot + 3._dl*F2**3*(F3dot**2 - 2._dl*F3*F3dotdot) + &
        4._dl*F2*F3**2*(-(F1dot*F3dot) + F3*F1dotdot + 3._dl*(-2._dl*F2dot + rho_d)*(F2dot + rho_d)) + &
        12._dl*F2**2*F3*(F3dot*(F2dot + rho_d) + F3*(F2dotdot - rho_d_dot))))))/ &
        (a*(3._dl*F2**2 + F1*F3)**2*(3._dl*F2**2 + 2._dl*F1*F3)) + &
        (9._dl*F2**6*(1._dl + 1._dl/(3._dl*sqrt(F2**4/(3._dl*F2**2 + 2._dl*F1*F3)**2)))* &
        (6._dl*F1*F2*F3*(3._dl*F2**2 + F1*F3)*adotoa* &
        (F2*F3*F1dot + F1*(F2*F3dot + 2._dl*F3*(-F2dot + rho_d))) + &
        a*(-3._dl*F2**4*F3**2*F1dot**2 + &
        2._dl*F1*F2**2*F3**2*(-(F3*F1dot**2) + 3._dl*F2**2*F1dotdot) + &
        2._dl*F1**3*F3*(-2._dl*F2**2*F3dot**2 + F2*F3*(F2*F3dotdot - 2._dl*F3dot*(-2._dl*F2dot + rho_d)) - &
        2._dl*F3**2*(F2dot**2 - F2dot*rho_d + rho_d**2 + F2*(F2dotdot - rho_d_dot))) + &
        F1**2*F2*(4._dl*F3**3*F1dot*F2dot + F2**3*(-9._dl*F3dot**2 + 6._dl*F3*F3dotdot) + &
        2._dl*F2*F3**2*(-(F1dot*F3dot) + F3*F1dotdot + 6._dl*(F2dot - rho_d)*rho_d) + &
        12._dl*F2**2*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))))))/ &
        (a*(3._dl*F2**2 + F1*F3)**2*(3._dl*F2**2 + 2._dl*F1*F3)))/(16._dl*F1**2*F2**4*F3**2)

    if (.not. (abs(mu1) < huge(1._dl) .and. abs(mu2) < huge(1._dl))) return
    eft_cache%EFT_mu1 = 2._dl*real(mu1, dl)
    eft_cache%EFT_mu2 = 2._dl*real(mu2, dl)
    if (present(resolved)) resolved = .true.

    ! The original transcription of arXiv:1609.03599 App. C (double precision, dust-only background, EFTcPP too large
    ! by about 1/adotoa; unusable: |mu_1|/H^2 up to 1e17 on healthy models), replaced above:
!   real(dl) :: F1, F2, F3, F1dot, F2dot, rho_d, rho_d_dot, F1dotdot, F2dotdot, F3dot, F3dotdot
!   real(dl) :: EFTcP, EFTcPP
    !
!   !< Define useful quantities used for the mass stability
!   rho_d = (eft_cache%grhob_t + eft_cache%grhoc_t)/a**2._dl
!   rho_d_dot = -3._dl*eft_cache%adotoa/a*rho_d
!   EFTcP = eft_cache%EFTcdot/(a*eft_cache%adotoa) + 2._dl*eft_cache%EFTc/a
    !
!   !< Define some EFT functions needed by stability
!   EFTcPP = (eft_cache%EFTcdotdot/eft_cache%adotoa + (3._dl - eft_cache%Hdot/eft_cache%adotoa)*eft_cache%EFTcdot + &
!       2._dl*eft_cache%adotoa*eft_cache%EFTc)/(a**2*eft_cache%adotoa**2)
    !
!   !< general expressions for Fs
!   F1 = 2._dl*(1._dl + eft_cache%EFTOmegaV) + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V
    !
!   F2 = eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V + &
!       (eft_cache%adotoa*(2._dl*(1._dl + eft_cache%EFTOmegaV) + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))/ &
!       a + eft_cache%adotoa*eft_cache%EFTOmegaP
    !
!   F3 = (2._dl*eft_cache%EFTc)/a**2. + 4._dl*eft_par_cache%h0_Mpc**2._dl*eft_cache%EFTGamma1V - &
!       (6._dl*eft_par_cache%h0_Mpc*eft_cache%adotoa*eft_cache%EFTGamma2V)/a - &
!       (3._dl*eft_cache%adotoa**2.* &
!       (2._dl*(1._dl + eft_cache%EFTOmegaV) + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))/a**2. - &
!       (6._dl*eft_cache%adotoa**2.*eft_cache%EFTOmegaP)/a
    !
!   F1dot = eft_cache%adotoa*(2._dl*eft_cache%EFTOmegaP + 3._dl*eft_cache%EFTGamma3P + eft_cache%EFTGamma4P)
    !
!   F2dot = ((2._dl + 2._dl*eft_cache%EFTOmegaV + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + &
!       a*eft_cache%EFTOmegaP)*eft_cache%Hdot + &
!       a**2._dl*eft_par_cache%h0_Mpc**2._dl*eft_cache%adotoa*eft_cache%EFTGamma2P + &
!       eft_cache%adotoa**2._dl*(-2._dl - 2._dl*eft_cache%EFTOmegaV - 3._dl*eft_cache%EFTGamma3V - &
!       eft_cache%EFTGamma4V + 2._dl*a*eft_cache%EFTOmegaP + 3._dl*a*eft_cache%EFTGamma3P + a*eft_cache%EFTGamma4P + &
!       a**2._dl*eft_cache%EFTOmegaPP))/a**2._dl
    !
!   F3dot = -((4._dl*eft_cache%EFTc*eft_cache%adotoa + &
!       6._dl*a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V*eft_cache%Hdot - &
!       2._dl*eft_cache%adotoa*(a*EFTcP - &
!       3._dl*(2._dl + 2._dl*eft_cache%EFTOmegaV + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + &
!       2._dl*a*eft_cache%EFTOmegaP)*eft_cache%Hdot + 2._dl*a**3._dl*eft_par_cache%h0_Mpc**2*eft_cache%EFTGamma1P) + &
!       6._dl*a*eft_par_cache%h0_Mpc*eft_cache%adotoa**2._dl*(-eft_cache%EFTGamma2V + a*eft_cache%EFTGamma2P) - &
!       3._dl*eft_cache%adotoa**3._dl* &
!       (4._dl + 4._dl*eft_cache%EFTOmegaV + 6._dl*eft_cache%EFTGamma3V + 2._dl*eft_cache%EFTGamma4V - &
!       3._dl*a*eft_cache%EFTGamma3P - a*eft_cache%EFTGamma4P - 2._dl*a**2*eft_cache%EFTOmegaPP))/a**3._dl)
    !
!   F1dotdot = (2._dl*eft_cache%EFTOmegaP*eft_cache%Hdot + &
!       eft_cache%Hdot*(3._dl*eft_cache%EFTGamma3P + eft_cache%EFTGamma4P) + &
!       a*eft_cache%adotoa**2._dl*(2._dl*eft_cache%EFTOmegaPP + 3._dl*eft_cache%EFTGamma3PP + eft_cache%EFTGamma4PP))/a
    !
!   F2dotdot = (a**2._dl*eft_par_cache%h0_Mpc**2._dl*eft_cache%Hdot*eft_cache%EFTGamma2P + &
!       eft_cache%adotoa*eft_cache%Hdot* &
!       (-8._dl - 8._dl*eft_cache%EFTOmegaV - 12._dl*eft_cache%EFTGamma3V - 4._dl*eft_cache%EFTGamma4V + &
!       5._dl*a*eft_cache%EFTOmegaP + 9._dl*a*eft_cache%EFTGamma3P + 3._dl*a*eft_cache%EFTGamma4P + &
!       3._dl*a**2._dl*eft_cache%EFTOmegaPP) + &
!       (2._dl + 2._dl*eft_cache%EFTOmegaV + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + &
!       a*eft_cache%EFTOmegaP)*eft_cache%Hdotdot + &
!       a**3._dl*eft_par_cache%h0_Mpc**2._dl*eft_cache%adotoa**2._dl*eft_cache%EFTGamma2PP + &
!       eft_cache%adotoa**3._dl*(4._dl + 4._dl*eft_cache%EFTOmegaV + 6._dl*eft_cache%EFTGamma3V + &
!       2._dl*eft_cache%EFTGamma4V - 4._dl*a*eft_cache%EFTOmegaP - 6._dl*a*eft_cache%EFTGamma3P - &
!       2._dl*a*eft_cache%EFTGamma4P + 2._dl*a**2._dl*eft_cache%EFTOmegaPP + 3._dl*a**2._dl*eft_cache%EFTGamma3PP + &
!       a**2._dl*eft_cache%EFTGamma4PP + a**3._dl*eft_cache%EFTOmegaPPP))/a**3._dl
    !
!   F3dotdot = (4._dl*eft_cache%EFTc*(3._dl*eft_cache%adotoa**2._dl - eft_cache%Hdot) + &
!       2._dl*(a*EFTcP*eft_cache%Hdot - &
!       3._dl*(2._dl + 2._dl*eft_cache%EFTOmegaV + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + &
!       2*a*eft_cache%EFTOmegaP)*eft_cache%Hdot**2 + &
!       2._dl*a**3._dl*eft_par_cache%h0_Mpc**2._dl*eft_cache%Hdot*eft_cache%EFTGamma1P - &
!       3._dl*a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V*eft_cache%Hdotdot) - &
!       6._dl*eft_cache%adotoa*(a*eft_par_cache%h0_Mpc*eft_cache%Hdot* &
!       (-4._dl*eft_cache%EFTGamma2V + 3._dl*a*eft_cache%EFTGamma2P) + &
!       (2._dl + 2._dl*eft_cache%EFTOmegaV + 3._dl*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + &
!       2._dl*a*eft_cache%EFTOmegaP)*eft_cache%Hdotdot) + &
!       eft_cache%adotoa**2._dl*(-8._dl*a*EFTcP + &
!       3._dl*eft_cache%Hdot*(24._dl + 24._dl*eft_cache%EFTOmegaV + 36._dl*eft_cache%EFTGamma3V + &
!       12._dl*eft_cache%EFTGamma4V + 4._dl*a*eft_cache%EFTOmegaP - 15._dl*a*eft_cache%EFTGamma3P - &
!       5._dl*a*eft_cache%EFTGamma4P - 10._dl*a**2._dl*eft_cache%EFTOmegaPP) + &
!       2._dl*a**2._dl*(EFTcPP + 2._dl*a**2._dl*eft_par_cache%h0_Mpc**2._dl*eft_cache%EFTGamma1PP)) - &
!       6._dl*a*eft_par_cache%h0_Mpc*eft_cache%adotoa**3* &
!       (2*eft_cache%EFTGamma2V + a*(-2._dl*eft_cache%EFTGamma2P + a*eft_cache%EFTGamma2PP)) - &
!       3._dl*eft_cache%adotoa**4._dl* &
!       (12._dl + 12._dl*eft_cache%EFTOmegaV + 18._dl*eft_cache%EFTGamma3V + 6._dl*eft_cache%EFTGamma4V - &
!       4._dl*a*eft_cache%EFTOmegaP - 12._dl*a*eft_cache%EFTGamma3P - 4._dl*a*eft_cache%EFTGamma4P - &
!       2._dl*a**2._dl*eft_cache%EFTOmegaPP + 3._dl*a**2._dl*eft_cache%EFTGamma3PP + a**2._dl*eft_cache%EFTGamma4PP + &
!       2._dl*a**3._dl*eft_cache%EFTOmegaPPP))/a**4._dl
    !
!   !< compute the mass eigenvalues
!   eft_cache%EFT_mu1 = ((-4._dl*F1*F2**2._dl*F3*(-2._dl*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2._dl)/ &
!       (3._dl*F2**2._dl + F1*F3) + &
!       (4._dl*F1*F3*(-2._dl*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2._dl)/(3._dl + (F1*F3)/F2**2._dl) - &
!       (9._dl*F2**6*(1._dl + 1._dl/(3.*sqrt(F2**4._dl/(3._dl*F2**2._dl + 2._dl*F1*F3)**2._dl)))* &
!       (6._dl*F1*F2*F3*(3._dl*F2**2._dl + F1*F3)*eft_cache%adotoa* &
!       (F2*F3*F1dot + F1*(F2*F3dot + 2._dl*F3*(-F2dot + rho_d))) + &
!       a*(-3._dl*F2**4._dl*F3**2._dl*F1dot**2._dl + &
!       2._dl*F1*F2**2._dl*F3**2._dl*(-(F3*F1dot**2._dl) + 3._dl*F2**2._dl*F1dotdot) + &
!       2._dl*F1**3._dl*F3*(-2._dl*F2**2._dl*F3dot**2._dl + F2*F3*(F2*F3dotdot - 2._dl*F3dot*(-2._dl*F2dot + rho_d)) - &
!       2._dl*F3**2._dl*(F2dot**2._dl - F2dot*rho_d + rho_d**2._dl + F2*(F2dotdot - rho_d_dot))) + &
!       F1**2._dl*F2*(4._dl*F3**3._dl*F1dot*F2dot + F2**3._dl*(-9._dl*F3dot**2._dl + 6._dl*F3*F3dotdot) + &
!       2._dl*F2*F3**2._dl*(-(F1dot*F3dot) + F3*F1dotdot + 6._dl*(F2dot - rho_d)*rho_d) + &
!       12._dl*F2**2._dl*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))))))/ &
!       (a*(3._dl*F2**2._dl + F1*F3)**2._dl*(3._dl*F2**2._dl + 2._dl*F1*F3)) + &
!       (2._dl*F1*F2**2._dl*(6._dl*F1*F2**2._dl*F3*(3._dl*F2**2._dl + F1*F3)*eft_cache%adotoa* &
!       (-2._dl*F3**2._dl*F1dot + 3._dl*F2**2._dl*F3dot + 6._dl*F2*F3*(-F2dot + rho_d)) + &
!       a*(6._dl*F2**4._dl*F3**3._dl*F1dot**2._dl + 2._dl*F1**3._dl*F3**2._dl*(-2._dl*F3*F2dot + F2*F3dot)**2._dl + &
!       F1**2._dl*F2*F3*(-8._dl*F3**3._dl*F1dot*F2dot + F2**3._dl*(-9._dl*F3dot**2._dl + 6._dl*F3*F3dotdot) + &
!       4._dl*F2*F3**2._dl*(F1dot*F3dot - F3*F1dotdot + 3._dl*(F2dot - rho_d)*rho_d) + &
!       12._dl*F2**2._dl*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))) + &
!       F1*(4._dl*F2**2._dl*F3**4._dl*F1dot**2._dl + 9._dl*F2**6._dl*(-3._dl*F3dot**2._dl + 2._dl*F3*F3dotdot) - &
!       12._dl*F2**4._dl*F3**2._dl*(F3*F1dotdot + 3._dl*rho_d*(-F2dot + rho_d)) + &
!       36._dl*F2**5._dl*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))))))/ &
!       (a*(3._dl*F2**2._dl + F1*F3)*(3._dl*F2**2._dl + 2._dl*F1*F3)))/(16._dl*F1**2._dl*F2**4._dl*F3**2._dl)
    !
!   eft_cache%EFT_mu2 = ((-4._dl*F1*F2**2._dl*F3*(-2._dl*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2._dl)/ &
!       (3._dl*F2**2._dl + F1*F3) + &
!       (4._dl*F1*F3*(-2._dl*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2._dl)/(3._dl + (F1*F3)/F2**2._dl) - &
!       (2._dl*F2**2._dl*F3*(6._dl*F1*F2**2._dl*F3*(3._dl*F2**2._dl + F1*F3)*eft_cache%adotoa* &
!       (9._dl*F2**4._dl*F1dot + 6._dl*F1*F2**2._dl*F3*F1dot + &
!       F1**2._dl*(-3._dl*F2**2._dl*F3dot + 2._dl*F3*(F3*F1dot + 3._dl*F2*(F2dot - rho_d)))) + &
!       a*(-27._dl*F2**8._dl*F3*F1dot**2._dl - 2._dl*F1**5._dl*F3**2._dl*(-2._dl*F3*F2dot + F2*F3dot)**2._dl + &
!       18._dl*F1*F2**6._dl*F3*(-2._dl*F3*F1dot**2._dl + 3._dl*F2**2._dl*F1dotdot) + &
!       18._dl*F1**2._dl*F2**4._dl*F3* &
!       (-(F1dot*(F3**2._dl*F1dot - 2._dl*F2*F3*F2dot + F2**2._dl*F3dot)) + 3._dl*F2**2._dl*F3*F1dotdot) + &
!       2._dl*F1**3._dl*F2**2._dl*(-2._dl*F3**4._dl*F1dot**2._dl + 12._dl*F2*F3**3._dl*F1dot*F2dot + &
!       9._dl*F2**4._dl*(F3dot**2._dl - F3*F3dotdot) + &
!       6._dl*F2**2._dl*F3**2._dl*(-3._dl*F2dot**2._dl - F1dot*F3dot + 2._dl*F3*F1dotdot - 3._dl*F2dot*rho_d + &
!       3._dl*rho_d**2._dl) + 18._dl*F2**3._dl*F3*(F3dot*rho_d + F3*(F2dotdot - rho_d_dot))) + &
!       F1**4._dl*F2*F3*(8._dl*F3**3._dl*F1dot*F2dot + 3._dl*F2**3._dl*(F3dot**2._dl - 2._dl*F3*F3dotdot) + &
!       4._dl*F2*F3**2._dl*(-(F1dot*F3dot) + F3*F1dotdot + 3._dl*(-2._dl*F2dot + rho_d)*(F2dot + rho_d)) + &
!       12._dl*F2**2._dl*F3*(F3dot*(F2dot + rho_d) + F3*(F2dotdot - rho_d_dot))))))/ &
!       (a*(3._dl*F2**2._dl + F1*F3)**2._dl*(3._dl*F2**2._dl + 2._dl*F1*F3)) + &
!       (9._dl*F2**6._dl*(1._dl + 1._dl/(3._dl*sqrt(F2**4._dl/(3._dl*F2**2._dl + 2._dl*F1*F3)**2._dl)))* &
!       (6._dl*F1*F2*F3*(3._dl*F2**2._dl + F1*F3)*eft_cache%adotoa* &
!       (F2*F3*F1dot + F1*(F2*F3dot + 2._dl*F3*(-F2dot + rho_d))) + &
!       a*(-3._dl*F2**4._dl*F3**2._dl*F1dot**2._dl + &
!       2._dl*F1*F2**2._dl*F3**2._dl*(-(F3*F1dot**2._dl) + 3._dl*F2**2._dl*F1dotdot) + &
!       2._dl*F1**3._dl*F3*(-2._dl*F2**2._dl*F3dot**2._dl + F2*F3*(F2*F3dotdot - 2._dl*F3dot*(-2._dl*F2dot + rho_d)) - &
!       2._dl*F3**2*(F2dot**2._dl - F2dot*rho_d + rho_d**2._dl + F2*(F2dotdot - rho_d_dot))) + &
!       F1**2*F2*(4._dl*F3**3._dl*F1dot*F2dot + F2**3._dl*(-9._dl*F3dot**2._dl + 6._dl*F3*F3dotdot) + &
!       2._dl*F2*F3**2._dl*(-(F1dot*F3dot) + F3*F1dotdot + 6._dl*(F2dot - rho_d)*rho_d) + &
!       12._dl*F2**2._dl*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))))))/ &
!       (a*(3._dl*F2**2._dl + F1*F3)**2._dl*(3._dl*F2**2._dl + 2._dl*F1*F3)))/(16._dl*F1**2._dl*F2**4._dl*F3**2._dl)

    end subroutine EFTCAMBModelComputeMassFactors

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the quasi-static expression of mu, from the general quasi-static form
    !! (1606.05339 App. A, the same as QS_mu_gen). The numerator and the denominator vanish together when pi
    !! does not couple to the metric: in that case the GR value mu = 1 is returned instead of 0/0.
    subroutine EFTCAMBModelComputeQSMu(self, a, k, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    real(dl), intent(in)                          :: k             !< the input scale.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: num, den, den_scale, slip_num, slip_den, slip_den_scale

    ! general quasi-static form (1606.05339 App. A). The former closed form had a wrong mass term (designer
    ! f(R): 4(Hddot - 3H^3) instead of Hddot - 2H^3) and no gamma3, gamma4, gamma6.
    call EFTCAMBQSGeneralTerms(a, k, eft_par_cache, eft_cache, num, den, den_scale, &
        slip_num, slip_den, slip_den_scale)
    if (EFTCAMBQSDegenerate(den, den_scale)) then
        eft_cache%QS_mu = 1._dl
    else
        eft_cache%QS_mu = num/den
    end if

    end subroutine EFTCAMBModelComputeQSMu

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the quasi-static expression of Sigma (lensing response).
    !! Sigma = mu (1 + gamma)/2 with mu the quasi-static mu and gamma = Phi/Psi the quasi-static
    !! gravitational slip, both from the general quasi-static form used by EFTCAMBModelComputeQSMuSigma
    !! (the same as QS_sigma_gen). A degenerate ratio (0/0, no coupling of pi to the metric) is set to its GR
    !! value.
    subroutine EFTCAMBModelComputeQSSigma(self, a, k, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    real(dl), intent(in)                          :: k             !< the input scale.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    real(dl) :: mu, gamma
    real(dl) :: mu_num, mu_den, mu_den_scale, slip_num, slip_den, slip_den_scale

    ! quasi-static mu and slip (general form):
    call EFTCAMBQSGeneralTerms(a, k, eft_par_cache, eft_cache, mu_num, mu_den, mu_den_scale, &
        slip_num, slip_den, slip_den_scale)
    if (EFTCAMBQSDegenerate(mu_den, mu_den_scale)) then
        mu = 1._dl
    else
        mu = mu_num/mu_den
    end if
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
    !> Subroutine that tests the model's own scalar ghost and gradient conditions. The generic conditions of the
    !! stability check extend up to GLPV: beyond it the check uses these instead (when EFT_ghost_stability or
    !! EFT_gradient_stability is on), and refuses the models that have none. The default has none.
    subroutine EFTCAMBModelBeyondGLPVStability(self, a, eft_par_cache, eft_cache, ghost, gradient, implemented, stable)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache
    logical, intent(in)                           :: ghost         !< whether to test the scalar ghost condition.
    logical, intent(in)                           :: gradient      !< whether to test the scalar gradient condition.
    logical, intent(out)                          :: implemented   !< whether the model has its own conditions.
    logical, intent(out)                          :: stable        !< whether the tested conditions are met.

    implemented = .false.
    stable = .true.

    end subroutine EFTCAMBModelBeyondGLPVStability

    ! ---------------------------------------------------------------------------------------------
    !> Function that reads model_background_a_final, the scale factor where a model background table ends, and
    !! returns its log. The table can extend past today: the models impose and read their conditions today at a = 1
    !! inside it. a_final < 1 is refused (today would be outside the table).
    function EFTCAMBModelReadBackgroundXFinal(self, Ini, a_final_default, eft_error) result(x_final)

    implicit none

    class(EFTCAMB_model)  :: self            !< the base class
    class(TIniFile)       :: Ini             !< input ini file
    real(dl), intent(in)  :: a_final_default !< the model's default a_final
    integer               :: eft_error       !< error code: 0 all fine, 1 initialization failed
    real(dl)              :: x_final         !< log(a_final)

    real(dl) :: a_final

    a_final = Ini%Read_Double('model_background_a_final', a_final_default)
    x_final = 0._dl
    if (.not. (a_final >= 1._dl)) then
        call GlobalError('EFTCAMB error: model_background_a_final must be >= 1 (today, a = 1, must be inside the ' // &
            'model background table)', error_unsupported_params)
        eft_error = 1
        return
    end if
    x_final = log(a_final)

    end function EFTCAMBModelReadBackgroundXFinal

    ! ---------------------------------------------------------------------------------------------
    !> Function that tests the phenomenological prior w_DE <= -1/3, enforced when EFT_wDE_prior = T.
    !! Models with a w_DE parametrization override it; the others have nothing to test.
    function EFTCAMBModelwDEPrior(self, a, eft_par_cache, eft_cache)

    implicit none

    class(EFTCAMB_model)                          :: self          !< the base class
    real(dl), intent(in)                          :: a             !< the input scale factor.
    !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
    !< the EFTCAMB timestep cache that contains all the physical values.
    type(TEFTCAMB_timestep_cache), intent(inout) :: eft_cache

    logical :: EFTCAMBModelwDEPrior !< true if the prior is met (or the model has no w_DE parametrization).

    EFTCAMBModelwDEPrior = .true.

    end function EFTCAMBModelwDEPrior

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the scale factor where the background solution of the model starts. Earlier, a model
    !! whose background is a table returns a frozen state: the perturbations should not be switched on there.
    !! The default 0 means that the background is defined at every a.
    function EFTCAMBModelBackgroundStart(self)

    implicit none

    class(EFTCAMB_model), intent(in) :: self      !< the base class
    real(dl) :: EFTCAMBModelBackgroundStart       !< the scale factor where the background solution starts.

    EFTCAMBModelBackgroundStart = 0._dl

    end function EFTCAMBModelBackgroundStart

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that checks that the evolution of pi is regular after the perturbation turn-on. A model whose pi
    !! variable is singular somewhere (e.g. Horndeski on the pi path: pi = delta_phi/phi' where phi' changes sign)
    !! returns the reason in message, and the run is refused; '' otherwise. Default: regular.
    subroutine EFTCAMBModelPiEvolutionCheck(self, turn_on, message)

    implicit none

    class(EFTCAMB_model), intent(in) :: self      !< the base class
    real(dl), intent(in)             :: turn_on   !< the scale factor of the perturbation turn-on.
    character(len=*), intent(out)    :: message   !< why the pi evolution is singular ('' if it is regular).

    message = ''

    end subroutine EFTCAMBModelPiEvolutionCheck

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

    ! (EFTCAMB_use_background is rejected for the other models when the parameters are set)
    EFTCAMBModelComputeDeltaTau = 0._dl
    call GlobalError('EFTCAMB error: EFTCAMB_use_background = True, EFTCAMBModelComputeDeltaTau only ' // &
        'supported by the Horndeski module (EFTflag=5)', error_unsupported_params)

    end function EFTCAMBModelComputeDeltaTau

    function EFTCAMBModelComputeDeltaCosmicT(self, a1, a2)

    implicit none

    class(EFTCAMB_model) :: self          !< the base class
    real(dl), intent(in) :: a1            !< the input scale factor.
    real(dl), intent(in) :: a2            !< the input scale factor.

    real(dl)             :: EFTCAMBModelComputeDeltaCosmicT !< the output t

    ! (EFTCAMB_use_background is rejected for the other models when the parameters are set)
    EFTCAMBModelComputeDeltaCosmicT = 0._dl
    call GlobalError('EFTCAMB error: EFTCAMB_use_background = True, EFTCAMBModelComputeDeltaCosmicT only ' // &
        'supported by the Horndeski module (EFTflag=5)', error_unsupported_params)

    end function EFTCAMBModelComputeDeltaCosmicT

    function EFTCAMBModelComputeSoundHorizon(self, a)

    implicit none

    class(EFTCAMB_model) :: self          !< the base class
    real(dl), intent(in) :: a             !< the input scale factor.

    real(dl)             :: EFTCAMBModelComputeSoundHorizon !< the output r_s

    ! (EFTCAMB_use_background is rejected for the other models when the parameters are set)
    EFTCAMBModelComputeSoundHorizon = 0._dl
    call GlobalError('EFTCAMB error: EFTCAMB_use_background = True, EFTCAMBModelComputeSoundHorizon only ' // &
        'supported by the Horndeski module (EFTflag=5)', error_unsupported_params)

    end function EFTCAMBModelComputeSoundHorizon

    end module EFTCAMB_abstract_model

! ----------------------------------------------------------------------------------------

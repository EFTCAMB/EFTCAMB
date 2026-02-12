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

!> @file 09_EFTCAMB_main.f90
!! This file contains the general EFTCAMB driver. We have a type that encapsulate all
!! EFTCAMB parameters and stuff and this is owned by CAMB to do all EFT calculations.


!----------------------------------------------------------------------------------------
!> This module contains the general EFTCAMB driver. We have a type that encapsulate all
!! EFTCAMB parameters and stuff and this is owned by CAMB to do all EFT calculations.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_main

    use precision
    use IniObjects
    use MpiUtils
    use classes
    use EFT_def
    use EFTCAMB_abstract_model
    use EFTCAMB_abstract_model_full
    use EFTCAMB_abstract_model_designer
    use EFTCAMB_pure_EFT_std
    use EFTCAMB_Reparametrized_Horndeski
    use EFTCAMB_Horndeski
    use EFTCAMB_designer_fR
    use EFTCAMB_designer_full3p1_wDE
    use EFTCAMB_designer_mc_quintessence
    use EFTCAMB_LE_Horava
    use EFTCAMB_Kmouflage_Mod
    use EFTCAMB_ADE_Mod
    use EFTCAMB_FM_quintessence
    use EFTCAMB_Omega_Lambda_gamma
    use EFTCAMB_Omega_Lambda_alpha
    use EFTCAMB_full_Beyond_Horndeski
    use EFTCAMB_full_Extended_Galileon
    use EFTCAMB_designer_ShiftSym_alphaB
    use EFTCAMB_full_Scaling_Cubic

    implicit none

    private

    public TEFTCAMB

    !----------------------------------------------------------------------------------------
    !> This is the main object for EFTCAMB and contains all the necessary ingredients to
    !! perform the EFT calculations. For the physical details of the calculations implemented
    !! please refer to the Numerical Notes (https://arxiv.org/abs/1405.3590).

    type, extends(TCambComponent) :: TEFTCAMB

        ! EFTCAMB model selection flags:
        integer   :: EFTflag              !< Main EFTCAMB model selection flag. Decides one of the four modes to run EFTCAMB.
        integer   :: PureEFTmodel         !< Model selection flag for pure EFT models.
        integer   :: AltParEFTmodel       !< Model selection flag for alternative EFT parametrizations.
        integer   :: DesignerEFTmodel     !< Model selection flag for designer mapping EFT models.
        integer   :: FullMappingEFTmodel  !< Model selection flag for full mapping EFT models.

        ! EFTCAMB stability flags:
        logical   :: EFT_ghost_math_stability      !< Flag that decides whether to enforce the mathematical ghost condition. This is less conservative than the physical requirement. Enforces the condition on both scalar and tensor sectors.
        logical   :: EFT_mass_math_stability       !< Flag that decides whether to enforce the mathematical mass condition. This should prevent fast growing instabilities but is not fully rigorous. Use is deprecated.
        !
        logical   :: EFT_ghost_stability           !< Flag that decides whether to enforce the physical ghost condition. Works up to Horndeski.
        logical   :: EFT_gradient_stability        !< Flag that decides whether to enforce the physical gradient condition. Works up to Horndeski.
        logical   :: EFT_mass_stability            !< Flag that decides whether to enforce the physical mass condition. Works up to Horndeski.
        real(dl)  :: EFT_mass_stability_rate       !< Flag that sets the rate for the mass instability in units of Hubble time.
        !
        logical   :: EFT_additional_priors         !< Flag that extablishes whether to use additional priors that are related to the specific model. Each model has the possibility of implementing their own additional stability requirements.
        
        ! Positvity bounds flags
        logical   :: EFT_positivity_bounds         !< Flag that decides whether to enforce the positivity bounds. Works up to Horndeski and for EFTFLAG=1 and (Omega,Lambda)-models.
        logical   :: EFT_minkowski_limit           !< Flag that checks whether there exists a healthy Minkowski limit. Works up to Horndeski.

        ! QSA flags
        logical   :: EFT_check_QSA             !< Flag that decides whether check QSA of mu
        real(dl)  :: EFT_QSA_threshold         !< Tolerence for relative difference between (mu,sigma) and that of QSA
        real(dl)  :: EFT_QSA_time              !< Minimum scale factor at which the code checks for QSA
        real(dl)  :: EFT_QSA_k                 !< Minimum k at which the code checks for QSA 

        ! Initial condition type
        integer   :: EFT_IC_type               !< Flag selecting initial condition types. 1: GR, 2: assume Omega=const. at IC time, 3: assume Omega=const. Gamma=const. at IC time

        ! EFTCAMB working flags:
        integer   :: EFTCAMB_feedback_level        !< Amount of feedback that is printed to screen.
        real(dl)  :: EFTCAMB_back_turn_on          !< Smallest scale factor that the code should
                                                   !! consider when computing the background. Set to zero (default), change background at all times.
        real(dl)  :: EFTCAMB_pert_turn_on          !< Smallest scale factor at which EFTCAMB evolves perturbations. Default set to a=0.01.
        real(dl)  :: EFTCAMB_GR_threshold          !< Tollerance on the dimensionless EFT functions to be considered GR. Default is 10^-8.
        real(dl)  :: EFTCAMB_stability_time        !< Minimum scale factor at which the code checks for stability of a theory
        real(dl)  :: EFTCAMB_stability_threshold   !< Threshold for the stability module to consider the model stable.
        logical   :: EFTCAMB_model_is_designer     !< Logical flag that establishes whether the model is designer or not.
        logical   :: EFTCAMB_effective_w0wa        !< Logical flag that establishes whether the model does have an effective w0wa parametrization.
        logical   :: EFTCAMB_use_background        !< Whether use the background conformal time, cosmic time and sound horizon computed by EFTCAMB instead of integrating them in CAMB. Currently only supported by the Horndeski module.
        logical   :: EFTCAMB_evolve_delta_phi      !< Flag that decides whether integrate delta_phi instead of pi for covariant theories.
        logical   :: EFTCAMB_evolve_metric_h       !< Flag that decides whether integrate h' instead of using constraints.
        logical   :: EFTCAMB_skip_stability        !< Flag that decides whether skip all stability checks.
        logical   :: EFTCAMB_skip_RGR              !< Flag that decides whether skip the return to GR checks.

        ! EFTCAMB output root:
        character(LEN=:), allocatable :: outroot   !< The root for auxiliary EFTCAMB output.

        ! EFTCAMB model:
        class(EFTCAMB_model), allocatable :: model !< This is the EFTCAMB model in the main class.

    contains

        ! utility functions:
        procedure :: EFTCAMB_init_from_file        => read_EFTCAMB_flags_from_file  !< subroutine that initializes EFTCAMB from an INI file.
        procedure :: EFTCAMB_init_model_from_file  => init_EFTCAMB_model_from_file  !< subroutine that initializes the selected model from file.
        procedure :: EFTCAMB_print_header          => print_EFTCAMB_header          !< subroutine that prints to screen the EFTCAMB header.
        procedure :: EFTCAMB_print_CosmoMC_header  => print_EFTCosmoMC_header       !< subroutine that prints to screen the EFTCosmoMC header.
        procedure :: EFTCAMB_print_model_feedback  => print_EFTCAMB_flags           !< subroutine that prints to screen the model flags and parameters.
        ! model allocation:
        procedure :: EFTCAMB_allocate_model            => allocate_EFTCAMB_model            !< subroutine that, based on the model selection flags allocates the EFTCAMB model.
        procedure :: EFTCAMB_read_model_selection      => read_EFTCAMB_model_selection      !< subroutine that reads the model selection parameters. Just a wrapper to the model specific subroutine.
        procedure :: EFTCAMB_allocate_model_functions  => allocate_EFTCAMB_model_functions  !< subroutine that, based on the model specific selection flags allocates the EFTCAMB model functions.
        procedure :: EFTCAMB_read_model_parameters     => read_EFTCAMB_model_parameters     !< subroutine that reads the model parameters. Just a wrapper to the model specific subroutine.
        procedure :: Effective_w_wa                    => EFTCAMB_effective_w0wa_values     !< subroutine that sets the effective w0wa values.
        ! python interface:
        procedure, nopass :: SelfPointer => TEFTCAMB_SelfPointer

    end type TEFTCAMB

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads from file the values of the EFTCAMB flags.
    subroutine read_EFTCAMB_flags_from_file( self, Ini )

        implicit none

        class(TEFTCAMB)  :: self  !< the base class
        type(TIniFile)   :: Ini   !< Input ini file

        character(len=:), allocatable :: outroot

        ! read from the INI file the main EFT flag:
        self%EFTflag              = Ini%Read_Int( 'EFTflag'            , 0 )
        if ( self%EFTflag == 0 ) return
        ! read the model selection flags:
        self%PureEFTmodel = 0
        self%AltParEFTmodel = 0
        self%DesignerEFTmodel = 0
        self%FullMappingEFTmodel = 0
        if ( self%EFTflag == 1 ) then
          self%PureEFTmodel         = Ini%Read_Int( 'PureEFTmodel'       , 1 )
        else if ( self%EFTflag == 2 ) then
          self%AltParEFTmodel       = Ini%Read_Int( 'AltParEFTmodel'     , 1 )
        else if ( self%EFTflag == 3 ) then
          self%DesignerEFTmodel     = Ini%Read_Int( 'DesignerEFTmodel'   , 1 )
        else if ( self%EFTflag == 4 ) then
          self%FullMappingEFTmodel  = Ini%Read_Int( 'FullMappingEFTmodel', 1 )
        end if

        ! mathematical stability:
        self%EFT_ghost_math_stability     = Ini%Read_Logical( 'EFT_ghost_math_stability', .false. )
        self%EFT_mass_math_stability      = Ini%Read_Logical( 'EFT_mass_math_stability' , .false. )

        ! physical stability:
        self%EFT_ghost_stability        = Ini%Read_Logical( 'EFT_ghost_stability'      , .true.  )
        self%EFT_gradient_stability     = Ini%Read_Logical( 'EFT_gradient_stability'   , .true.  )
        self%EFT_mass_stability         = Ini%Read_Logical( 'EFT_mass_stability'       , .false. )
        self%EFT_mass_stability_rate    = Ini%Read_Double ( 'EFT_mass_stability_rate'  , 10._dl  )

        ! model specific priors:
        self%EFT_additional_priors      = Ini%Read_Logical( 'EFT_additional_priors'    , .true. )

        !) positivity bounds:
        self%EFT_positivity_bounds      = Ini%Read_Logical( 'EFT_positivity_bounds'    , .false. )

        !) Existence of a healthy Minkowski limit:
        self%EFT_minkowski_limit        = Ini%Read_Logical( 'EFT_minkowski_limit'    , .false. )

        !) initial condition types
        self%EFT_IC_type                = Ini%Read_Int( 'EFT_IC_type', 1 )
        
        ! EFTCAMB working flags:
        self%EFTCAMB_feedback_level      = Ini%Read_Int    ( 'feedback_level'             , 1      )
        self%EFTCAMB_back_turn_on        = Ini%Read_Double ( 'EFTCAMB_back_turn_on'       , 1.d-8  )
        self%EFTCAMB_pert_turn_on        = Ini%Read_Double ( 'EFTCAMB_turn_on_time'       , 1.d-2  )
        self%EFTCAMB_GR_threshold        = Ini%Read_Double ( 'EFTCAMB_GR_threshold'       , 1.d-8  )
        self%EFTCAMB_stability_time      = Ini%Read_Double ( 'EFTCAMB_stability_time'     , 1.d-10 )
        self%EFTCAMB_stability_threshold = Ini%Read_Double ( 'EFTCAMB_stability_threshold', 0._dl  )
        self%EFTCAMB_effective_w0wa      = Ini%Read_Logical( 'EFTCAMB_effective_w0wa'     , .false.)
        ! flags currently only supported by the Horndeski module
        self%EFTCAMB_use_background      = Ini%Read_Logical( 'EFTCAMB_use_background'     , .false. )
        if ( self%EFTCAMB_use_background .and. ( self%EFTflag /= 5 ) )then
            write(*,*) 'EFTCAMB Error: EFTCAMB_use_background=True currently only supported by the Horndeski module (EFTflag=5)'
            call MpiStop("EFTCAMB error")
        end if
        self%EFTCAMB_evolve_delta_phi         = Ini%Read_Logical( 'EFTCAMB_evolve_delta_phi', .false. )
        if ( self%EFTCAMB_evolve_delta_phi .and. ( self%EFTflag /= 5) ) then
            write(*,*) "EFTCAMB Error: EFTCAMB_evolve_delta_phi=True currently only supported by the Horndeski module (EFTflag = 5)."
            call MpiStop('EFTCAMB error')
        end if
        self%EFTCAMB_evolve_metric_h          = Ini%Read_Logical( 'EFTCAMB_evolve_metric_h', .false. )
        self%EFTCAMB_skip_stability           = Ini%Read_Logical('EFTCAMB_skip_stability', .false.)
        self%EFTCAMB_skip_RGR                 = Ini%Read_Logical('EFTCAMB_skip_RGR', .false.)

        ! Output root for debug purposes:
        outroot = Ini%Read_String('output_root')
        if (outroot == '') outroot = 'EFTCAMB'
        outroot = trim(outroot) // '_'
        self%outroot = trim( outroot )

    end subroutine read_EFTCAMB_flags_from_file

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the selected model from file.
    subroutine init_EFTCAMB_model_from_file( self, Ini, eft_error )

        implicit none

        class(TEFTCAMB) :: self      !< the base class
        type(TIniFile)  :: Ini       !< input ini file
        integer         :: eft_error !< error code: 0 all fine, 1 initialization failed

        eft_error = 0
        ! allocate model:
        call self%EFTCAMB_allocate_model(eft_error)
        if ( eft_error > 0 ) return
        ! read the parameters defining the model from file:
        call self%EFTCAMB_read_model_selection( Ini, eft_error )
        if ( eft_error > 0 ) return
        ! allocate model functions and parameters:
        call self%EFTCAMB_allocate_model_functions( Ini, eft_error )
        if ( eft_error > 0 ) return
        ! read model parameters from file:
        call self%EFTCAMB_read_model_parameters( Ini, eft_error )
        if ( eft_error > 0 ) return
        ! compute model number of parameters:
        call self%model%compute_param_number()

    end subroutine init_EFTCAMB_model_from_file

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that gets the effective w0wa values.
    subroutine EFTCAMB_effective_w0wa_values( self, w0, wa )

        implicit none

        class(TEFTCAMB) :: self      !< the base class
        real(dl)        :: w0        !< the effective w0 value (output)
        real(dl)        :: wa        !< the effective wa value (output)

        if ( self%EFTCAMB_effective_w0wa ) then
            w0 = self%model%w0
            wa = self%model%wa
        else
            w0 = -1._dl
            wa = 0._dl
        end if

    end subroutine EFTCAMB_effective_w0wa_values


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the EFTCAMB header.
    subroutine print_EFTCAMB_header( self )

        implicit none

        class(TEFTCAMB) :: self       !< the base class

        ! check feedback level:
        if ( .not. self%EFTCAMB_feedback_level > 0 ) return
        ! if GR return:
        if ( self%EFTflag == 0 ) return
        ! print the header:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') "     ______________________   __  ______  "
        write(*,'(a)') "    / __/ __/_  __/ ___/ _ | /  |/  / _ ) "
        write(*,'(a)') "   / _// _/  / / / /__/ __ |/ /|_/ / _  | "
        write(*,'(a)') "  /___/_/   /_/  \___/_/ |_/_/  /_/____/  "//" "//EFTCAMB_version
        write(*,'(a)') "  "
        write(*,'(a)') "***************************************************************"

    end subroutine print_EFTCAMB_header

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the EFTCosmoMC header.
    subroutine print_EFTCosmoMC_header( self )

        implicit none

        class(TEFTCAMB) :: self       !< the base class

        ! check feedback level:
        if ( .not. self%EFTCAMB_feedback_level > 0 ) return
        ! if GR return:
        if ( self%EFTflag == 0 ) return
        ! print the header:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') "     ___________________                    __  ________"
        write(*,'(a)') "    / __/ __/_  __/ ___/__  ___ __ _  ___  /  |/  / ___/"
        write(*,'(a)') "   / _// _/  / / / /__/ _ \(_-</  ' \/ _ \/ /|_/ / /__  "
        write(*,'(a)') "  /___/_/   /_/  \___/\___/___/_/_/_/\___/_/  /_/\___/  "
        write(*,'(a)') "  "
        write(*,'(a)') "  "//EFTCAMB_version
        write(*,'(a)') "***************************************************************"

    end subroutine print_EFTCosmoMC_header

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen informations about the model and the model parameters.
    subroutine print_EFTCAMB_flags( self, print_params )

        implicit none

        class(TEFTCAMB)     :: self         !< the base class
        logical, optional   :: print_params !< optional flag that decised whether to print numerical values
                                            !! of the parameters.

        character(len=500)  :: temp_name
        real(dl)            :: temp_value
        integer             :: i

        ! if GR return:
        if ( self%EFTflag == 0 ) return

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            write(*,*) 'EFTCAMB WARNING: trying to call EFTCAMB model feedback without allocating the model'
            call MpiStop('EFTCAMB error')
        end if

        ! check feedback level:
        if ( .not. self%EFTCAMB_feedback_level > 0 ) return
        ! print settings flags:
        write(*,*)
        write(*,*) 'EFTCAMB settings:'
        write(*,*) ' EFTCAMB_feedback_level      =', self%EFTCAMB_feedback_level
        write(*,*) ' EFTCAMB_back_turn_on        =', self%EFTCAMB_back_turn_on
        write(*,*) ' EFTCAMB_pert_turn_on        =', self%EFTCAMB_pert_turn_on
        write(*,*) ' EFTCAMB_GR_threshold        =', self%EFTCAMB_GR_threshold
        write(*,*) ' EFTCAMB_stability_time      =', self%EFTCAMB_stability_time
        write(*,*) ' EFTCAMB_stability_threshold =', self%EFTCAMB_stability_threshold
        write(*,*) ' EFTCAMB_effective_w0wa      =', self%EFTCAMB_effective_w0wa
        
        ! print stability flags:
        write(*,*)
        write(*,*) 'EFTCAMB stability flags:'
        write(*,*) ' Math ghost stability   = ', self%EFT_ghost_math_stability
        write(*,*) ' Math mass stability    = ', self%EFT_mass_math_stability
        if ( self%EFT_mass_math_stability ) then
            write(*,*) ' The use of the math mass condition is deprecated.'
            write(*,*) ' Make sure you know what you are doing...'
        end if
        write(*,*) ' Ghost stability        = ', self%EFT_ghost_stability
        write(*,*) ' Gradient stability     = ', self%EFT_gradient_stability
        write(*,*) ' Mass stability         = ', self%EFT_mass_stability
        write(*,'(a,E12.5)') '  Mass stability rate    = ', self%EFT_mass_stability_rate
        write(*,*) ' Additional priors      = ', self%EFT_additional_priors
        write(*,*) ' Positivity bounds      = ', self%EFT_positivity_bounds
        write(*,*) ' Minkowski limit      = ', self%EFT_minkowski_limit
        if ( self%EFT_positivity_bounds ) then
            write(*,*) ' The exact calculation of the positivity bounds is implemented'
            write(*,*) ' only for models with an abstract parametrization using gamma functions:'
            write(*,*) "  PureEFTmodel (EFTFlag = 1) and OmegaLambda with Gamma's (EFTFlag=2, AltParEFTmodel=2)"  
        end if
        write(*,*)
        ! print model selection flags:
        write(*,*)              'EFTCAMB model flags:'
        write(*,"(A24,I3)")     '   EFTflag             =', self%EFTflag
        if ( self%EFTflag == 1 ) &
            write(*,"(A24,I3)") '   PureEFTmodel        =', self%PureEFTmodel
        if ( self%EFTflag == 2 ) &
            write(*,"(A24,I3)") '   AltParEFTmodel      =', self%AltParEFTmodel
        if ( self%EFTflag == 3 ) &
            write(*,"(A24,I3)") '   DesignerEFTmodel    =', self%DesignerEFTmodel
        if ( self%EFTflag == 4 ) &
            write(*,"(A24,I3)") '   FullMappingEFTmodel =', self%FullMappingEFTmodel
        ! print model informations:
        call self%model%feedback( print_params )
        ! leave one white line:
        write(*,*)

    end subroutine print_EFTCAMB_flags

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the EFTCAMB model based on the model selection flags.
    !! If implementing a new model this is the place to allocate it.
    subroutine allocate_EFTCAMB_model( self, eft_error )

        implicit none

        class(TEFTCAMB) :: self       !< the base class
        integer         :: eft_error  !< error code: 0 all fine, 1 initialization failed

        ! check the allocation of the model:
        if ( allocated(self%model) ) deallocate(self%model)

        ! do the allocation:
        select case ( self%EFTflag )

            case (0)     ! GR: no need to allocate

            case (1)     ! Pure EFT:

                select case ( self%PureEFTmodel )
                    case(1)
                        allocate( EFTCAMB_std_pure_EFT::self%model )
                        call self%model%init( 'Standard Pure EFT', 'Standard Pure EFT' )
                    case default
                        if (self%EFTCAMB_feedback_level > 0) then
                            write(*,'(a,I3)') 'No model corresponding to EFTFlag =', self%EFTflag
                            write(*,'(a,I3)') 'and PureEFTmodel =', self%PureEFTmodel
                            write(*,'(a)')    'Please select an appropriate model:'
                            write(*,'(a)')    'PureEFTmodel=1  standard Pure EFT'
                        end if
                        eft_error = 1
                        return
                end select

            case (2)     ! Alternative EFT:

                select case ( self%AltParEFTmodel )
                    case(1)
                        allocate( EFTCAMB_RPH::self%model )
                        call self%model%init( 'RPH', 'RPH' )
                    case(2) !Lambda_Omega with gamma functions
                        allocate( EFTCAMB_OmegaLambda_gamma::self%model )
                        call self%model%init( 'OL gamma', 'OL gamma' )
                    case(3) !Lambda_Omega with alpha functions
                        allocate( EFTCAMB_OmegaLambda_alpha::self%model )
                        call self%model%init( 'OL alpha', 'OL alpha' )
                    case(4) !Shift Simmetric Gravity with alpha_B
                        allocate( EFTCAMB_ShiftSym_alphaB::self%model )
                        call self%model%init( 'Shift Symmetric alpha B', 'Shift Symmetric alpha B' )
                    case default
                        if (self%EFTCAMB_feedback_level > 0) then
                            write(*,'(a,I3)') 'No model corresponding to EFTFlag =', self%EFTflag
                            write(*,'(a,I3)') 'and AltParEFTmodel =', self%AltParEFTmodel
                            write(*,'(a)')    'Please select an appropriate model:'
                            write(*,'(a)')    'AltParEFTmodel=1  reparametrized Horndeski'
                            write(*,'(a)')    'AltParEFTmodel=2  Lambda-Omega parametrization with gamma functions'
                            write(*,'(a)')    'AltParEFTmodel=3  Lambda-Omega parametrization with alpha functions'
                            write(*,'(a)')    'AltParEFTmodel=4  Shift Symmetric Gravity with alpha_B'
                        end if
                        eft_error = 1
                        return
                end select

            case (3)     ! Designer mapping EFT:

                select case ( self%DesignerEFTmodel )
                    case(1)
                        allocate( EFTCAMB_fR_designer::self%model )
                        call self%model%init( 'Designer f(R)', 'Designer f(R)' )
                    case(2)
                        allocate( EFTCAMB_des_mc_quint::self%model )
                        call self%model%init( 'Designer minimally coupled quintessence', 'Designer minimally coupled quintessence' )
                    case(3)
                        allocate( EFTCAMB_full3p1_wDE::self%model )
                        call self%model%init( 'Designer full3p1', 'Designer full3p1' )
                    case default
                        if (self%EFTCAMB_feedback_level > 0) then
                            write(*,'(a,I3)') 'No model corresponding to EFTFlag =', self%EFTflag
                            write(*,'(a,I3)') 'and DesignerEFTmodel =', self%DesignerEFTmodel
                            write(*,'(a)')    'Please select an appropriate model:'
                            write(*,'(a)')    'DesignerEFTmodel=1  designer f(R)'
                            write(*,'(a)')    'DesignerEFTmodel=2  designer minimally coupled quintessence'
                            write(*,'(a)')    'DesignerEFTmodel=3  designer full 3+1 theory'
                        end if
                        eft_error = 1
                        return
                end select

            case (4)     ! Full mapping EFT:

                select case ( self%FullMappingEFTmodel )
                    case(1)
                        allocate( EFTCAMB_Horava::self%model )
                        call self%model%init( 'Horava', 'Horava' )
                    case(2)
                        allocate( EFTCAMB_ADE::self%model )
                        call self%model%init( 'Acoustic Dark Energy', 'Acoustic Dark Energy' )
                    case(3)
                        allocate( EFTCAMB_Kmouflage::self%model )
                        call self%model%init( 'K-mouflage', 'K-mouflage' )
                    case(4)
                        allocate( EFTCAMB_5e::self%model )
                        call self%model%init( 'Quintessence', 'Quintessence' )
                    case(5)
                        allocate( EFTCAMB_Beyond_Horndeski::self%model)
                        call self%model%init( 'Beyond Horndeski', 'Beyond Horndeski' )
                    case(6)
                        allocate( EFTCAMB_Scaling_Cubic::self%model)
                        call self%model%init( 'Scaling Cubic Galileon', 'Scaling Cubic Galileon' )
                    case(7)
                        allocate( EFTCAMB_Extended_Galileon::self%model)
                        call self%model%init( 'Extended Galileon', 'Extended Galileon' )
                    case default
                        if (self%EFTCAMB_feedback_level > 0) then
                            write(*,'(a,I3)') 'No model corresponding to EFTFlag =', self%EFTflag
                            write(*,'(a,I3)') 'and FullMappingEFTmodel =', self%FullMappingEFTmodel
                            write(*,'(a)')    'Please select an appropriate model:'
                            write(*,'(a)')    'FullMappingEFTmodel=1  Horava gravity'
                            write(*,'(a)')    'FullMappingEFTmodel=2  Acoustic Dark Energy'
                            write(*,'(a)')    'FullMappingEFTmodel=3  K-mouflage'
                            write(*,'(a)')    'FullMappingEFTmodel=4  Quintessence'
                            write(*,'(a)')    'FullMappingEFTmodel=5  Beyond Horndeski'
                            write(*,'(a)')    'FullMappingEFTmodel=6  Scaling Cubic Galileon'
                            write(*,'(a)')    'FullMappingEFTmodel=7  Extended Galileon'
                        end if
                        eft_error = 1
                        return
                end select
            
            case (5)
                allocate( EFTCAMB_Hdsk::self%model )
                call self%model%init( 'Horndeski', 'Horndeski' )


            case default ! not found:
                if (self%EFTCAMB_feedback_level > 0) then
                    write(*,'(a,I3)') 'No model corresponding to EFTFlag =', self%EFTflag
                    write(*,'(a)') 'Please select an appropriate model:'
                    write(*,'(a)') 'EFTFlag=0  GR code'
                    write(*,'(a)') 'EFTFlag=1  Pure EFT'
                    write(*,'(a)') 'EFTFlag=2  EFT alternative parametrizations'
                    write(*,'(a)') 'EFTFlag=3  designer mapping EFT'
                    write(*,'(a)') 'EFTFlag=4  full mapping EFT'
                    write(*,'(a)') 'EFTFlag=5  Horndeski'
                end if
                eft_error = 1
                return

        end select

        ! now store the designer flag:
        select type ( model => self%model )
            class is ( EFTCAMB_full_model )
            self%EFTCAMB_model_is_designer = .False.
            class is ( EFTCAMB_designer_model )
            self%EFTCAMB_model_is_designer = .True.
        end select

    end subroutine allocate_EFTCAMB_model

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the model selection parameters. Just a wrapper to the model specific subroutine.
    subroutine read_EFTCAMB_model_selection( self, Ini, eft_error )

        implicit none

        class(TEFTCAMB) :: self      !< the base class
        type(TIniFile)  :: Ini       !< Input ini file
        integer         :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            if (self%EFTCAMB_feedback_level > 0) then
                write(*,*) 'EFTCAMB WARNING: trying to call EFTCAMB model read_model_selection'
                write(*,*) ' without allocating the model'
            end if
            eft_error = 1
            return
        end if

        ! call the model specific read parameters:
        call self%model%read_model_selection( Ini, eft_error )

    end subroutine read_EFTCAMB_model_selection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that, based on the model specific selection flags allocates the EFTCAMB model functions.
    !! Just a wrapper to the model specific subroutine.
    subroutine allocate_EFTCAMB_model_functions( self, Ini, eft_error )

        implicit none

        class(TEFTCAMB) :: self       !< the base class
        type(TIniFile)  :: Ini        !< Input ini file
        integer         :: eft_error  !< error code: 0 all fine, 1 initialization failed

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            if (self%EFTCAMB_feedback_level > 0) then
                write(*,*) 'EFTCAMB WARNING: trying to call EFTCAMB model allocate_model_selection'
                write(*,*) ' without allocating the model'
            end if
            eft_error = 1
            return
        end if

        ! call the model specific read parameters:
        call self%model%allocate_model_selection( Ini, eft_error )

    end subroutine allocate_EFTCAMB_model_functions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the model parameters. Just a wrapper to the model specific subroutine.
    subroutine read_EFTCAMB_model_parameters( self, Ini, eft_error )

        implicit none

        class(TEFTCAMB) :: self       !< the base class
        type(TIniFile)  :: Ini        !< Input ini file
        integer         :: eft_error  !< error code: 0 all fine, 1 initialization failed

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            if (self%EFTCAMB_feedback_level > 0) then
                write(*,*) 'EFTCAMB WARNING: trying to call EFTCAMB model read_model_parameters_from_file'
                write(*,*) ' without allocating the model'
            end if
            eft_error = 1
            return
        end if

        ! call the model specific read parameters:
        call self%model%init_model_parameters_from_file( Ini, eft_error )

        ! now store the effective w0wa flag:
        if (self%EFTCAMB_effective_w0wa .and. self%model%effective_w0wa) then
            self%EFTCAMB_effective_w0wa = .True.
        else
            self%EFTCAMB_effective_w0wa = .False.
        end if

    end subroutine read_EFTCAMB_model_parameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine for the python interface:
    subroutine TEFTCAMB_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TEFTCAMB), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P
        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TEFTCAMB_SelfPointer

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_main

!----------------------------------------------------------------------------------------

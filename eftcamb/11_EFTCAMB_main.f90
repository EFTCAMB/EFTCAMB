!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2016 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 11_EFTCAMB_main.f90
!! This file contains the general EFTCAMB driver. We have a type that encapsulate all
!! EFTCAMB parameters and stuff and this is owned by CAMB to do all EFT calculations.


!----------------------------------------------------------------------------------------
!> This module contains the general EFTCAMB driver. We have a type that encapsulate all
!! EFTCAMB parameters and stuff and this is owned by CAMB to do all EFT calculations.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_main

    use precision
    use IniFile
    use AMLutils
    use EFT_def
    use EFTCAMB_abstract_model
    use EFTCAMB_pure_EFT_std
    use EFTCAMB_Reparametrized_Horndeski
    use EFTCAMB_designer_fR
    use EFTCAMB_designer_mc_quintessence
    use EFTCAMB_LE_Horava

    implicit none

    private

    public EFTCAMB_version, EFTCAMB

    character(LEN=*), parameter :: EFTCAMB_version = 'V3.0 Oct16'

    !----------------------------------------------------------------------------------------
    !> This is the main object for EFTCAMB. Get one of these and you can use all the stuff
    !! in EFTCAMB.
    type EFTCAMB

        ! EFTCAMB output root:
        character(LEN=:), allocatable :: outroot !< The root for auxiliary EFTCAMB output.

        ! EFTCAMB model selection flags:
        integer   :: EFTflag              !< Main EFTCAMB model selection flag. Decides one of the four modes to run EFTCAMB.
        integer   :: PureEFTmodel         !< Model selection flag for pure EFT models.
        integer   :: AltParEFTmodel       !< Model selection flag for alternative EFT parametrizations.
        integer   :: DesignerEFTmodel     !< Model selection flag for designer mapping EFT models.
        integer   :: FullMappingEFTmodel  !< Model selection flag for full mapping EFT models.

        ! EFTCAMB stability flags:
        logical   :: EFT_mathematical_stability  !< Flag that extablishes wether to use mathematical stability.
        logical   :: EFT_physical_stability      !< Flag that extablishes wether to use physical stability.
        logical   :: EFT_additional_priors       !< Flag that extablishes wether to use some additional priors that are related to the specific model.

        ! EFTCAMB model:
        class(EFTCAMB_model), allocatable :: model !< This is the EFTCAMB model in the main class.

        ! EFTCAMB working flags:
        integer   :: EFTCAMB_feedback_level      !< Amount of feedback that is printed to screen.
        real(dl)  :: EFTCAMB_turn_on_time        !< Scale factor at which EFTCAMB becomes active. Default set to EFTturnonpiInitial in 01_EFT_def.f90.

    contains

        ! utility functions:
        procedure :: EFTCAMB_init_from_file        => read_EFTCAMB_flags_from_file  !< subroutine that initializes EFTCAMB from an INI file.
        procedure :: EFTCAMB_init_model_from_file  => init_EFTCAMB_model_from_file  !< subroutine that initializes the selected model from file.
        procedure :: EFTCAMB_print_header          => print_EFTCAMB_header          !< subroutine that prints to screen the EFTCAMB header.
        procedure :: EFTCAMB_print_model_feedback  => print_EFTCAMB_flags           !< subroutine that prints to screen the model flags and parameters.
        ! model allocation:
        procedure :: EFTCAMB_allocate_model            => allocate_EFTCAMB_model            !< subroutine that, based on the model selection flags allocates the EFTCAMB model.
        procedure :: EFTCAMB_read_model_selection      => read_EFTCAMB_model_selection      !< subroutine that reads the model selection parameters. Just a wrapper to the model specific subroutine.
        procedure :: EFTCAMB_allocate_model_functions  => allocate_EFTCAMB_model_functions  !< subroutine that, based on the model specific selection flags allocates the EFTCAMB model functions.
        procedure :: EFTCAMB_read_model_parameters     => read_EFTCAMB_model_parameters     !< subroutine that reads the model parameters. Just a wrapper to the model specific subroutine.

    end type EFTCAMB

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads from file the values of the EFTCAMB flags.
    subroutine read_EFTCAMB_flags_from_file( self, Ini )

        implicit none

        class(EFTCAMB)      :: self       !< the base class
        type(TIniFile)      :: Ini        !< Input ini file

        ! read from the INI file the main EFT flag:
        self%EFTflag              = Ini_Read_Int_File( Ini, 'EFTflag'            , 0 )
        ! read the model selection flags:
        self%PureEFTmodel         = Ini_Read_Int_File( Ini, 'PureEFTmodel'       , 0 )
        self%AltParEFTmodel       = Ini_Read_Int_File( Ini, 'AltParEFTmodel'     , 0 )
        self%DesignerEFTmodel     = Ini_Read_Int_File( Ini, 'DesignerEFTmodel'   , 0 )
        self%FullMappingEFTmodel  = Ini_Read_Int_File( Ini, 'FullMappingEFTmodel', 0 )
        ! read the stability flags:
        self%EFT_mathematical_stability = Ini_Read_Logical_File( Ini, 'EFT_mathematical_stability', .true. )
        self%EFT_physical_stability     = Ini_Read_Logical_File( Ini, 'EFT_physical_stability'    , .true. )
        self%EFT_additional_priors      = Ini_Read_Logical_File( Ini, 'EFT_additional_priors'     , .true. )

        ! EFTCAMB working stuff:
        self%EFTCAMB_feedback_level     = Ini_Read_Int_File( Ini, 'feedback_level', 1 )
        self%EFTCAMB_turn_on_time       = Ini_Read_Double_File( Ini, 'EFTCAMB_turn_on_time', EFTturnonpiInitial )

    end subroutine read_EFTCAMB_flags_from_file

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the selected model from file.
    subroutine init_EFTCAMB_model_from_file( self, Ini )

        implicit none

        class(EFTCAMB)      :: self       !< the base class
        type(TIniFile)      :: Ini        !< Input ini file

        ! allocate model:
        call self%EFTCAMB_allocate_model()
        ! read the parameters defining the model from file:
        call self%EFTCAMB_read_model_selection( Ini )
        ! allocate model functions and parameters:
        call self%EFTCAMB_allocate_model_functions( )
        ! read model parameters from file:
        call self%EFTCAMB_read_model_parameters( Ini )
        ! compute model number of parameters:
        call self%model%compute_param_number()

    end subroutine init_EFTCAMB_model_from_file

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the EFTCAMB header.
    subroutine print_EFTCAMB_header( self )

        implicit none

        class(EFTCAMB) :: self       !< the base class

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
    !> Subroutine that prints to screen informations about the model and the model parameters.
    subroutine print_EFTCAMB_flags( self )

        implicit none

        class(EFTCAMB)      :: self       !< the base class

        character(len=500)  :: temp_name
        real(dl)            :: temp_value
        integer             :: i

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            write(*,*) 'EFTCAMB WARNING: trying to call EFTCAMB model feedback without allocating the model'
            call MpiStop('EFTCAMB error')
        end if
        ! check feedback level:
        if ( .not. self%EFTCAMB_feedback_level > 0 ) return
        ! if GR return:
        if ( self%EFTflag == 0 ) return
        ! print feedback flag:
        write(*,*)
        write(*,'(a, I3)') ' EFTCAMB feedback level  =', self%EFTCAMB_feedback_level

        ! print stability flags:
        write(*,*)
        write(*,*) 'EFTCAMB stability flags:'

        write(*,*) ' Mathematical stability = ', self%EFT_mathematical_stability
        write(*,*) ' Physical stability     = ', self%EFT_physical_stability
        write(*,*) ' Additional priors      = ', self%EFT_additional_priors
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
        call self%model%feedback( )
        ! leave one white line:
        write(*,*)

    end subroutine print_EFTCAMB_flags

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the EFTCAMB model based on the model selection flags.
    !! If implementing a new model this is the place to allocate it.
    subroutine allocate_EFTCAMB_model( self )

        implicit none

        class(EFTCAMB)      :: self       !< the base class

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
                        write(*,'(a,I3)') 'No model corresponding to EFTFlag =', self%EFTflag
                        write(*,'(a,I3)') 'and PureEFTmodel =', self%PureEFTmodel
                        write(*,'(a)')    'Please select an appropriate model:'
                        write(*,'(a)')    'PureEFTmodel=1  standard Pure EFT'
                        call MpiStop('EFTCAMB error')
                end select

            case (2)     ! Alternative EFT:

                select case ( self%AltParEFTmodel )
                    case(1)
                        allocate( EFTCAMB_RPH::self%model )
                        call self%model%init( 'RPH', 'RPH' )
                    case default
                        write(*,'(a,I3)') 'No model corresponding to EFTFlag =', self%EFTflag
                        write(*,'(a,I3)') 'and AltParEFTmodel =', self%AltParEFTmodel
                        write(*,'(a)')    'Please select an appropriate model:'
                        write(*,'(a)')    'AltParEFTmodel=1  reparametrized Horndeski'
                        call MpiStop('EFTCAMB error')
                end select

            case (3)     ! Designer mapping EFT:

                select case ( self%DesignerEFTmodel )
                    case(1)
                        allocate( EFTCAMB_fR_designer::self%model )
                        call self%model%init( 'Designer f(R)', 'Designer f(R)' )
                    case(2)
                        allocate( EFTCAMB_des_mc_quint::self%model )
                        call self%model%init( 'Designer minimally coupled quintessence', 'Designer minimally coupled quintessence' )
                    case default
                        write(*,'(a,I3)') 'No model corresponding to EFTFlag =', self%EFTflag
                        write(*,'(a,I3)') 'and DesignerEFTmodel =', self%DesignerEFTmodel
                        write(*,'(a)')    'Please select an appropriate model:'
                        write(*,'(a)')    'DesignerEFTmodel=1  designer f(R)'
                        write(*,'(a)')    'DesignerEFTmodel=2  designer minimally coupled quintessence'
                        call MpiStop('EFTCAMB error')
                end select

            case (4)     ! Full mapping EFT:

                select case ( self%FullMappingEFTmodel )
                    case(1)
                        allocate( EFTCAMB_Horava::self%model )
                        call self%model%init( 'Horava', 'Horava' )
                    case default
                        write(*,'(a,I3)') 'No model corresponding to EFTFlag =', self%EFTflag
                        write(*,'(a,I3)') 'and FullMappingEFTmodel =', self%FullMappingEFTmodel
                        write(*,'(a)')    'Please select an appropriate model:'
                        write(*,'(a)')    'FullMappingEFTmodel=1  Horava gravity'
                        call MpiStop('EFTCAMB error')
                end select

            case default ! not found:

                write(*,'(a,I3)') 'No model corresponding to EFTFlag =', self%EFTflag
                write(*,'(a)') 'Please select an appropriate model:'
                write(*,'(a)') 'EFTFlag=0  GR code'
                write(*,'(a)') 'EFTFlag=1  Pure EFT'
                write(*,'(a)') 'EFTFlag=2  EFT alternative parametrizations'
                write(*,'(a)') 'EFTFlag=3  designer mapping EFT'
                write(*,'(a)') 'EFTFlag=4  full mapping EFT'
                call MpiStop('EFTCAMB error')

        end select

    end subroutine allocate_EFTCAMB_model

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the model selection parameters. Just a wrapper to the model specific subroutine.
    subroutine read_EFTCAMB_model_selection( self, Ini )

        implicit none

        class(EFTCAMB)      :: self       !< the base class
        type(TIniFile)      :: Ini        !< Input ini file

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            write(*,*) 'EFTCAMB WARNING: trying to call EFTCAMB model read_model_selection'
            write(*,*) ' without allocating the model'
            call MpiStop('EFTCAMB error')
        end if

        ! call the model specific read parameters:
        call self%model%read_model_selection( Ini )

    end subroutine read_EFTCAMB_model_selection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that, based on the model specific selection flags allocates the EFTCAMB model functions.
    !! Just a wrapper to the model specific subroutine.
    subroutine allocate_EFTCAMB_model_functions( self )

        implicit none

        class(EFTCAMB)      :: self       !< the base class

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            write(*,*) 'EFTCAMB WARNING: trying to call EFTCAMB model allocate_model_selection'
            write(*,*) ' without allocating the model'
            call MpiStop('EFTCAMB error')
        end if

        ! call the model specific read parameters:
        call self%model%allocate_model_selection( )

    end subroutine allocate_EFTCAMB_model_functions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the model parameters. Just a wrapper to the model specific subroutine.
    subroutine read_EFTCAMB_model_parameters( self, Ini )

        implicit none

        class(EFTCAMB)      :: self       !< the base class
        type(TIniFile)      :: Ini        !< Input ini file

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            write(*,*) 'EFTCAMB WARNING: trying to call EFTCAMB model read_model_parameters_from_file'
            write(*,*) ' without allocating the model'
            call MpiStop('EFTCAMB error')
        end if

        ! call the model specific read parameters:
        call self%model%init_model_parameters_from_file( Ini )

    end subroutine read_EFTCAMB_model_parameters

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_main

!----------------------------------------------------------------------------------------

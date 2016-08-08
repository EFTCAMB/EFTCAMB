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

!> @file 9_EFTCAMB_main.f90
!! This file contains the general EFTCAMB driver. We have a type that encapsulate all
!! EFTCAMB parameters and stuff and this is owned by CAMB to do all EFT calculations.


!----------------------------------------------------------------------------------------
!> This module contains the general EFTCAMB driver. We have a type that encapsulate all
!! EFTCAMB parameters and stuff and this is owned by CAMB to do all EFT calculations.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_main

    use precision
    use EFTCAMB_abstract_model
    use IniFile

    implicit none

    character(LEN=*), parameter :: EFTCAMB_version = 'Ago16'

    !----------------------------------------------------------------------------------------
    !> This is the main object for EFTCAMB. Get one of these and you can use all the stuff
    !! in EFTCAMB.
    type EFTCAMB

        ! EFTCAMB model selection flags:
        integer   :: EFTflag              !< Main EFTCAMB model selection flag. Decides one of the four modes to run EFTCAMB.
        integer   :: PureEFTmodel         !< Model selection flag for pure EFT models.
        integer   :: AltParEFTmodel       !< Model selection flag for alternative EFT parametrizations.
        integer   :: DesignerEFTmodel     !< Model selection flag for designer mapping EFT models.
        integer   :: FullMappingEFTmodel  !< Model selection flag for full mapping EFT models.

        ! EFTCAMB stability flags:
        logical   :: EFT_mathematical_stability
        logical   :: EFT_physical_stability
        logical   :: EFTAdditionalPriors
        logical   :: MinkowskyPriors

        ! EFTCAMB model:
        type(EFTCAMB_model), allocatable :: model  !< This is the EFTCAMB model in the main class.

        ! EFTCAMB working flags:
        integer   :: EFTCAMB_feedback_level        !< Amount of feedback that is printed to screen.

    contains

        ! utility functions:
        procedure :: EFTCAMB_init_from_file       => read_EFTCAMB_flags_from_file  !< subroutine that initializes EFTCAMB from an INI file.
        procedure :: EFTCAMB_print_model_feedback => print_EFTCAMB_flags           !< subroutine that prints to screen the model flags and parameters.
        ! model allocation:
        procedure :: EFTCAMB_allocate_model => allocate_EFTCAMB_model        !< subroutine that, based on the model selection flags allocates the EFTCAMB model.

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
        self%EFTAdditionalPriors        = Ini_Read_Logical_File( Ini, 'EFTAdditionalPriors'       , .true. )
        self%MinkowskyPriors            = Ini_Read_Logical_File( Ini, 'MinkowskyPriors'           , .true. )

    end subroutine read_EFTCAMB_flags_from_file

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
            ! GR:
            case (0)
                allocate( self%model )
            ! Pure EFT:
            case (1)
            ! Alternative EFT:
            case (2)
            ! Designer mapping EFT:
            case (3)
            ! Full mapping EFT:
            case (4)
            ! not found:
            case default
                write(*,*) 'No model corresponding to EFTFlag=', self%EFTflag
                write(*,*) 'Please select an appropriate model:'
                write(*,*) 'EFTFlag=0  GR code'
                write(*,*) 'EFTFlag=1  Pure EFT'
                write(*,*) 'EFTFlag=2  EFT alternative parametrizations'
                write(*,*) 'EFTFlag=3  designer mapping EFT'
                write(*,*) 'EFTFlag=4  full mapping EFT'
                stop
        end select

    end subroutine allocate_EFTCAMB_model

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen informations about the model and the model parameters.
    subroutine print_EFTCAMB_flags( self )

        implicit none

        class(EFTCAMB)      :: self       !< the base class

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            write(*,*) 'EFTCAMB WARNING: trying to call EFTCAMB model feedback without allocationg the model'
            stop
        end if
        ! check feedback level:
        if ( .not. self%EFTCAMB_feedback_level > 0 ) return
        ! if GR return:
        if ( self%EFTflag == 0 ) return
        ! print stability flags:
        write(*,*)
        write(*,*) 'EFTCAMB stability flags:'

        write(*,*) ' Mathematical stability = ', self%EFT_mathematical_stability
        write(*,*) ' Physical stability     = ', self%EFT_physical_stability
        write(*,*) ' Additional priors      = ', self%EFTAdditionalPriors
        write(*,*) ' Minkowsky priors       = ', self%MinkowskyPriors
        write(*,*)
        ! print model selection flags:
        write(*,*)              'EFTCAMB model flags:'
        write(*,"(A24,I3)")     '   EFTflag             =', self%EFTflag
        if ( self%EFTflag == 1 ) &
            write(*,"(A24,I3)") '   PureEFTmodel        =', self%DesignerEFTmodel
        if ( self%EFTflag == 2 ) &
            write(*,"(A24,I3)") '   AltParEFTmodel      =', self%DesignerEFTmodel
        if ( self%EFTflag == 3 ) &
            write(*,"(A24,I3)") '   DesignerEFTmodel    =', self%DesignerEFTmodel
        if ( self%EFTflag == 4 ) &
            write(*,"(A24,I3)") '   FullMappingEFTmodel =', self%DesignerEFTmodel

    end subroutine print_EFTCAMB_flags

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_main

!----------------------------------------------------------------------------------------

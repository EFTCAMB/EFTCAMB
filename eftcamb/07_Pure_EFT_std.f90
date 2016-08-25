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

!> @file 07_Pure_EFT_std.f90
!! This file contains the definition of the Pure EFT model in which the EFT is described
!! by six functions of time and w_DE. Please refer to the numerical notes for details.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Pure EFT model in which the EFT is described
!! by six functions of time and w_DE. Please refer to the numerical notes for details.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_pure_EFT_std

    use precision
    use IniFile
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_neutral_parametrization_1D
    use EFTCAMB_constant_parametrization_1D
    use EFTCAMB_linear_parametrizations_1D
    use EFTCAMB_abstract_model_designer

    implicit none

    !----------------------------------------------------------------------------------------
    !> This is the pure EFT model with the six functions of time and w_DE.
    type, extends ( EFTCAMB_designer_model ) :: EFTCAMB_std_pure_EFT

        ! the pure EFT functions model selection flags:
        integer  :: PureEFTmodelOmega   !< Model selection flag for Pure EFT Omega.
        integer  :: EFTwDE              !< Model selection flag for Pure EFT w DE.
        integer  :: PureEFTmodelGamma1  !< Model selection flag for Pure EFT Gamma1.
        integer  :: PureEFTmodelGamma2  !< Model selection flag for Pure EFT Gamma2.
        integer  :: PureEFTmodelGamma3  !< Model selection flag for Pure EFT Gamma3.
        integer  :: PureEFTmodelGamma4  !< Model selection flag for Pure EFT Gamma4.
        integer  :: PureEFTmodelGamma5  !< Model selection flag for Pure EFT Gamma5.
        integer  :: PureEFTmodelGamma6  !< Model selection flag for Pure EFT Gamma6.

        ! selection flag for Horndeski:
        logical  :: PureEFTHorndeski    !< Selects wether to use the Horndeski bound on EFT functions.

        ! the pure EFT functions:
        class( parametrized_function_1D ), allocatable :: PureEFTOmega    !< The pure EFT function Omega.
        class( parametrized_function_1D ), allocatable :: PureEFTwDE      !< The pure EFT function w_DE.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma1   !< The pure EFT function Gamma1.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma2   !< The pure EFT function Gamma2.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma3   !< The pure EFT function Gamma3.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma4   !< The pure EFT function Gamma4.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma5   !< The pure EFT function Gamma5.
        class( parametrized_function_1D ), allocatable :: PureEFTGamma6   !< The pure EFT function Gamma6.

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBPureEFTstdReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBPureEFTstdAllocateModelSelection      !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBPureEFTstdInitModelParameters         !< subroutine taht initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBPureEFTstdInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.
        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBPureEFTstdComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBPureEFTstdFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBPureEFTstdParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBPureEFTstdParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBPureEFTstdParameterValues            !< subroutine that returns the i-th parameter value.
        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBPureEFTstdBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBPureEFTstdSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_dtauda                    => EFTCAMBPureEFTstdComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure :: compute_adotoa                    => EFTCAMBPureEFTstdComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBPureEFTstdComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.

    end type EFTCAMB_std_pure_EFT

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBPureEFTstdReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        ! read model selection flags:
        self%PureEFTmodelOmega  = Ini_Read_Int_File( Ini, 'PureEFTmodelOmega'  , 0 )
        self%EFTwDE             = Ini_Read_Int_File( Ini, 'EFTwDE'             , 0 )
        self%PureEFTmodelGamma1 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma1' , 0 )
        self%PureEFTmodelGamma2 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma2' , 0 )
        self%PureEFTmodelGamma3 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma3' , 0 )
        self%PureEFTmodelGamma4 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma4' , 0 )
        self%PureEFTmodelGamma5 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma5' , 0 )
        self%PureEFTmodelGamma6 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma6' , 0 )
        ! read the Horndeski flag:
        self%PureEFTHorndeski   = Ini_Read_Logical_File( Ini, 'PureEFTHorndeski' , .false. )

    end subroutine EFTCAMBPureEFTstdReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBPureEFTstdAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class

        ! allocate Omega:
        if ( allocated(self%PureEFTOmega) ) deallocate(self%PureEFTOmega)
        select case ( self%PureEFTmodelOmega )
            case(0)
                allocate( zero_parametrization_1D::self%PureEFTOmega )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTOmega )
            case(2)
                allocate( linear_parametrization_1D::self%PureEFTOmega )
            case default
                write(*,'(a,I3)') 'No model corresponding to PureEFTmodelOmega =', self%PureEFTmodelOmega
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate wDE:
        if ( allocated(self%PureEFTwDE) ) deallocate(self%PureEFTwDE)
        select case ( self%EFTwDE )
            case(0)
                allocate( wDE_LCDM_parametrization_1D::self%PureEFTwDE )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTwDE )
            case default
                write(*,'(a,I3)') 'No model corresponding to EFTwDE =', self%EFTwDE
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate Gamma1:
        if ( allocated(self%PureEFTGamma1) ) deallocate(self%PureEFTGamma1)
        select case ( self%PureEFTmodelGamma1 )
            case(0)
                allocate( zero_parametrization_1D::self%PureEFTGamma1 )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTGamma1 )
            case default
                write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma1 =', self%PureEFTmodelGamma1
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate Gamma2:
        if ( allocated(self%PureEFTGamma2) ) deallocate(self%PureEFTGamma2)
        select case ( self%PureEFTmodelGamma2 )
            case(0)
                allocate( zero_parametrization_1D::self%PureEFTGamma2 )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTGamma2 )
            case default
                write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma2 =', self%PureEFTmodelGamma2
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate Gamma3:
        if ( allocated(self%PureEFTGamma3) ) deallocate(self%PureEFTGamma3)
        select case ( self%PureEFTmodelGamma3 )
            case(0)
                allocate( zero_parametrization_1D::self%PureEFTGamma3 )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTGamma3 )
            case default
                write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma3 =', self%PureEFTmodelGamma3
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate the other functions only if not Horndeski:
        if ( .not. self%PureEFTHorndeski ) then
            ! allocate Gamma4:
            if ( allocated(self%PureEFTGamma4) ) deallocate(self%PureEFTGamma4)
            select case ( self%PureEFTmodelGamma4 )
                case(0)
                    allocate( zero_parametrization_1D::self%PureEFTGamma4 )
                case(1)
                    allocate( constant_parametrization_1D::self%PureEFTGamma4 )
                case default
                    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma4 =', self%PureEFTmodelGamma4
                    write(*,'(a)')    'Please select an appropriate model.'
            end select
            ! allocate Gamma5:
            if ( allocated(self%PureEFTGamma5) ) deallocate(self%PureEFTGamma5)
            select case ( self%PureEFTmodelGamma5 )
                case(0)
                    allocate( zero_parametrization_1D::self%PureEFTGamma5 )
                case(1)
                    allocate( constant_parametrization_1D::self%PureEFTGamma5 )
                case default
                    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma5 =', self%PureEFTmodelGamma5
                    write(*,'(a)')    'Please select an appropriate model.'
            end select
            ! allocate Gamma6:
            if ( allocated(self%PureEFTGamma6) ) deallocate(self%PureEFTGamma6)
            select case ( self%PureEFTmodelGamma6 )
                case(0)
                    allocate( zero_parametrization_1D::self%PureEFTGamma6 )
                case(1)
                    allocate( constant_parametrization_1D::self%PureEFTGamma6 )
                case default
                    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma6 =', self%PureEFTmodelGamma6
                    write(*,'(a)')    'Please select an appropriate model.'
            end select
        end if

        ! initialize the names:
        call self%PureEFTOmega%init ( 'EFTOmega' , '\Omega'       )
        call self%PureEFTwDE%init   ( 'EFTw'     , 'w'            )
        call self%PureEFTGamma1%init( 'EFTGamma1', '\gamma^{(1)}' )
        call self%PureEFTGamma2%init( 'EFTGamma2', '\gamma^{(2)}' )
        call self%PureEFTGamma3%init( 'EFTGamma3', '\gamma^{(3)}' )

        if ( .not. self%PureEFTHorndeski ) then
            call self%PureEFTGamma4%init( 'EFTGamma4', '\gamma^{(4)}' )
            call self%PureEFTGamma5%init( 'EFTGamma5', '\gamma^{(5)}' )
            call self%PureEFTGamma6%init( 'EFTGamma6', '\gamma^{(6)}' )
        end if

    end subroutine EFTCAMBPureEFTstdAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine taht initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBPureEFTstdInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.

    end subroutine EFTCAMBPureEFTstdInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBPureEFTstdInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        call self%PureEFTOmega%init_from_file ( Ini )
        call self%PureEFTwDE%init_from_file   ( Ini )
        call self%PureEFTGamma1%init_from_file( Ini )
        call self%PureEFTGamma2%init_from_file( Ini )
        call self%PureEFTGamma3%init_from_file( Ini )

        if ( .not. self%PureEFTHorndeski ) then
            call self%PureEFTGamma4%init_from_file( Ini )
            call self%PureEFTGamma5%init_from_file( Ini )
            call self%PureEFTGamma6%init_from_file( Ini )
        end if

    end subroutine EFTCAMBPureEFTstdInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBPureEFTstdComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_std_pure_EFT)  :: self   !< the base class

        self%parameter_number = 0
        self%parameter_number = self%parameter_number +self%PureEFTOmega%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTwDE%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTGamma1%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTGamma2%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTGamma3%parameter_number
        if ( .not. self%PureEFTHorndeski ) then
            self%parameter_number = self%parameter_number +self%PureEFTGamma4%parameter_number
            self%parameter_number = self%parameter_number +self%PureEFTGamma5%parameter_number
            self%parameter_number = self%parameter_number +self%PureEFTGamma6%parameter_number
        end if

    end subroutine EFTCAMBPureEFTstdComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBPureEFTstdFeedback( self )

        implicit none

        class(EFTCAMB_std_pure_EFT)  :: self   !< the base class

        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        if ( self%PureEFTHorndeski ) then
            write(*,"(a)")  '   Pure EFT Horndeski'
        end if
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
        ! print model functions informations:
        write(*,*)
        if ( self%PureEFTmodelOmega  /= 0 ) write(*,'(a,I3)') '   PureEFTmodelOmega   =', self%PureEFTmodelOmega
        if ( self%EFTwDE             /= 0 ) write(*,'(a,I3)') '   EFTwDE              =', self%EFTwDE
        if ( self%PureEFTmodelGamma1 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma1  =', self%PureEFTmodelGamma1
        if ( self%PureEFTmodelGamma2 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma2  =', self%PureEFTmodelGamma2
        if ( self%PureEFTmodelGamma3 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma3  =', self%PureEFTmodelGamma3
        if ( .not. self%PureEFTHorndeski ) then
            if ( self%PureEFTmodelGamma4 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma4  =', self%PureEFTmodelGamma4
            if ( self%PureEFTmodelGamma5 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma5  =', self%PureEFTmodelGamma5
            if ( self%PureEFTmodelGamma6 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma6  =', self%PureEFTmodelGamma6
        end if

        write(*,*)

        call self%PureEFTOmega%feedback()
        call self%PureEFTwDE%feedback()
        call self%PureEFTGamma1%feedback()
        call self%PureEFTGamma2%feedback()
        call self%PureEFTGamma3%feedback()
        if ( .not. self%PureEFTHorndeski ) then
            call self%PureEFTGamma4%feedback()
            call self%PureEFTGamma5%feedback()
            call self%PureEFTGamma6%feedback()
        end if

    end subroutine EFTCAMBPureEFTstdFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        stop 'IW: EFTCAMBPureEFTstdParameterNames'

    end subroutine EFTCAMBPureEFTstdParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

        stop 'IW: EFTCAMBPureEFTstdParameterNamesLatex'

    end subroutine EFTCAMBPureEFTstdParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

        stop 'IW: EFTCAMBPureEFTstdParameterValues'

    end subroutine EFTCAMBPureEFTstdParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBPureEFTstdBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%EFTOmegaV    = self%PureEFTOmega%value(a)
        eft_cache%EFTOmegaP    = self%PureEFTOmega%first_derivative(a)
        eft_cache%EFTOmegaPP   = self%PureEFTOmega%second_derivative(a)
        eft_cache%EFTOmegaPPP  = self%PureEFTOmega%third_derivative(a)
        eft_cache%EFTc         = ( eft_cache%adotoa**2 - eft_cache%Hdot )*( eft_cache%EFTOmegaV + 0.5_dl*a*eft_cache%EFTOmegaP ) &
            & -0.5_dl*( a*eft_cache%adotoa )**2*eft_cache%EFTOmegaPP&
            & +0.5_dl*eft_cache%grhov_t*( 1._dl+ self%PureEFTwDE%value(a) )
        eft_cache%EFTLambda    = -eft_cache%EFTOmegaV*( 2._dl*eft_cache%Hdot +eft_cache%adotoa**2 ) &
            & -a*eft_cache%EFTOmegaP*( 2._dl*eft_cache%adotoa**2 + eft_cache%Hdot ) &
            & -( a*eft_cache%adotoa )**2*eft_cache%EFTOmegaPP &
            & +self%PureEFTwDE%value(a)*eft_cache%grhov_t
        eft_cache%EFTcdot      = +0.5_dl*eft_cache%adotoa*eft_cache%grhov_t*( -3._dl*(1._dl +self%PureEFTwDE%value(a))**2 + a*self%PureEFTwDE%first_derivative(a) ) &
            & -eft_cache%EFTOmegaV*( eft_cache%Hdotdot -4._dl*eft_cache%adotoa*eft_cache%Hdot +2._dl*eft_cache%adotoa**3 ) &
            & +0.5_dl*a*eft_cache%EFTOmegaP*( -eft_cache%Hdotdot +eft_cache%adotoa*eft_cache%Hdot +eft_cache%adotoa**3) &
            & +0.5_dl*a**2*eft_cache%adotoa*eft_cache%EFTOmegaPP*( eft_cache%adotoa**2 -3._dl*eft_cache%Hdot ) &
            & -0.5_dl*(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP
        eft_cache%EFTLambdadot = -2._dl*eft_cache%EFTOmegaV*( eft_cache%Hdotdot -eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3 ) &
            & -a*eft_cache%EFTOmegaP*( 4._dl*eft_cache%adotoa*eft_cache%Hdot +eft_cache%Hdotdot ) &
            & -a**2*eft_cache%EFTOmegaPP*eft_cache%adotoa*( +2._dl*eft_cache%adotoa**2 +3._dl*eft_cache%Hdot )&
            & -(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP &
            & +eft_cache%grhov_t*eft_cache%adotoa*( a*self%PureEFTwDE%first_derivative(a) -3._dl*self%PureEFTwDE%value(a)*(1._dl +self%PureEFTwDE%value(a) ))

    end subroutine EFTCAMBPureEFTstdBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBPureEFTstdSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%EFTGamma1V  = self%PureEFTGamma1%value(a)
        eft_cache%EFTGamma1P  = self%PureEFTGamma1%first_derivative(a)
        eft_cache%EFTGamma2V  = self%PureEFTGamma2%value(a)
        eft_cache%EFTGamma2P  = self%PureEFTGamma2%first_derivative(a)
        eft_cache%EFTGamma3V  = self%PureEFTGamma3%value(a)
        eft_cache%EFTGamma3P  = self%PureEFTGamma3%first_derivative(a)
        if ( self%PureEFTHorndeski ) then
            eft_cache%EFTGamma4V  = -self%PureEFTGamma3%value(a)
            eft_cache%EFTGamma4P  = -self%PureEFTGamma3%first_derivative(a)
            eft_cache%EFTGamma4PP = -self%PureEFTGamma3%second_derivative(a)
            eft_cache%EFTGamma5V  = +0.5_dl*self%PureEFTGamma3%value(a)
            eft_cache%EFTGamma5P  = +0.5_dl*self%PureEFTGamma3%first_derivative(a)
            eft_cache%EFTGamma6V  = 0._dl
            eft_cache%EFTGamma6P  = 0._dl
        else
            eft_cache%EFTGamma4V  = self%PureEFTGamma4%value(a)
            eft_cache%EFTGamma4P  = self%PureEFTGamma4%first_derivative(a)
            eft_cache%EFTGamma4PP = self%PureEFTGamma4%second_derivative(a)
            eft_cache%EFTGamma5V  = self%PureEFTGamma5%value(a)
            eft_cache%EFTGamma5P  = self%PureEFTGamma5%first_derivative(a)
            eft_cache%EFTGamma6V  = self%PureEFTGamma6%value(a)
            eft_cache%EFTGamma6P  = self%PureEFTGamma6%first_derivative(a)
        end if

    end subroutine EFTCAMBPureEFTstdSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBPureEFTstdComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBPureEFTstdComputeDtauda               !< the output dtauda

        real(dl) :: temp

        temp = eft_cache%grhoa2 +eft_par_cache%grhov*a*a*self%PureEFTwDE%integral(a)
        EFTCAMBPureEFTstdComputeDtauda = sqrt(3/temp)

    end function EFTCAMBPureEFTstdComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBPureEFTstdComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%grhov_t = eft_par_cache%grhov*self%PureEFTwDE%integral(a)
        eft_cache%adotoa  = sqrt( ( eft_cache%grhom_t +eft_cache%grhov_t )/3._dl )

    end subroutine EFTCAMBPureEFTstdComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBPureEFTstdComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_std_pure_EFT)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%gpiv_t  = self%PureEFTwDE%value(a)*eft_cache%grhov_t
        eft_cache%Hdot    = -0.5_dl*( eft_cache%adotoa**2 +eft_cache%gpresm_t +eft_cache%gpiv_t )

        ! IW POSSIBLE BUG
        !eft_cache%Hdotdot = eft_cache%adotoa*( ( eft_cache%grhob_t +eft_cache%grhoc_t)/6._dl +2._dl*( eft_cache%grhor_t +eft_cache%grhog_t)/3._dl ) &
        !    & +eft_cache%adotoa*eft_cache%grhov_t*( 1._dl/6._dl +self%PureEFTwDE%value(a) +1.5_dl*self%PureEFTwDE%value(a)**2 -0.5_dl*a*self%PureEFTwDE%first_derivative(a) ) &
        !    & +eft_cache%adotoa*eft_cache%grhonu_tot/6._dl -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot

        eft_cache%Hdotdot = 2._dl*eft_cache%adotoa*eft_cache%Hdot &
                & + 0.5_dl*eft_cache%adotoa*(eft_cache%grhob_t + eft_cache%grhoc_t + 8._dl*(eft_cache%grhog_t+eft_cache%grhor_t)/3._dl)&
                & + 0.5_dl*eft_cache%adotoa*eft_cache%grhov_t*( (1._dl+self%PureEFTwDE%value(a) )*(1._dl+3._dl*self%PureEFTwDE%value(a)) -a*self%PureEFTwDE%first_derivative(a))&
                & + eft_cache%adotoa/6._dl*eft_cache%grhonu_tot -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot

    end subroutine EFTCAMBPureEFTstdComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_pure_EFT_std

!----------------------------------------------------------------------------------------

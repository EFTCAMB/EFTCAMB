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

!> @file 06_Pure_EFT_std.f90
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
    use EFTCAMB_constant_parametrization_1D
    use EFTCAMB_linear_parametrizations_1D
    use EFTCAMB_abstract_model

    implicit none

    !----------------------------------------------------------------------------------------
    !> This is the pure EFT model with the six functions of time and w_DE.
    type, extends ( EFTCAMB_model ) :: EFTCAMB_std_pure_EFT

        ! the pure EFT functions model selection flags:
        integer  :: PureEFTmodelOmega   !< Model selection flag for Pure EFT Omega.
        integer  :: EFTwDE              !< Model selection flag for Pure EFT w DE.
        integer  :: PureEFTmodelGamma1  !< Model selection flag for Pure EFT Gamma1.
        integer  :: PureEFTmodelGamma2  !< Model selection flag for Pure EFT Gamma2.
        integer  :: PureEFTmodelGamma3  !< Model selection flag for Pure EFT Gamma3.
        integer  :: PureEFTmodelGamma4  !< Model selection flag for Pure EFT Gamma4.
        integer  :: PureEFTmodelGamma5  !< Model selection flag for Pure EFT Gamma5.
        integer  :: PureEFTmodelGamma6  !< Model selection flag for Pure EFT Gamma6.

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
        procedure :: compute_dtauda                    => EFTCAMBPureEFTstdComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2)

    end type EFTCAMB_std_pure_EFT

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBPureEFTstdReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        self%PureEFTmodelOmega  = Ini_Read_Int_File( Ini, 'PureEFTmodelOmega'  , 0 )
        self%EFTwDE             = Ini_Read_Int_File( Ini, 'EFTwDE'             , 0 )
        self%PureEFTmodelGamma1 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma1' , 0 )
        self%PureEFTmodelGamma2 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma2' , 0 )
        self%PureEFTmodelGamma3 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma3' , 0 )
        self%PureEFTmodelGamma4 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma4' , 0 )
        self%PureEFTmodelGamma5 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma5' , 0 )
        self%PureEFTmodelGamma6 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma6' , 0 )

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
                allocate( parametrized_function_1D::self%PureEFTOmega )
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
                allocate( parametrized_function_1D::self%PureEFTwDE )
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
                allocate( parametrized_function_1D::self%PureEFTGamma1 )
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
                allocate( parametrized_function_1D::self%PureEFTGamma2 )
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
                allocate( parametrized_function_1D::self%PureEFTGamma3 )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTGamma3 )
            case default
                write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma3 =', self%PureEFTmodelGamma3
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate Gamma4:
        if ( allocated(self%PureEFTGamma4) ) deallocate(self%PureEFTGamma4)
        select case ( self%PureEFTmodelGamma4 )
            case(0)
                allocate( parametrized_function_1D::self%PureEFTGamma4 )
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
                allocate( parametrized_function_1D::self%PureEFTGamma5 )
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
                allocate( parametrized_function_1D::self%PureEFTGamma6 )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTGamma6 )
            case default
                write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma6 =', self%PureEFTmodelGamma6
                write(*,'(a)')    'Please select an appropriate model.'
        end select

        ! initialize the names:
        call self%PureEFTOmega%init ( 'EFTOmega' , '\Omega'       )
        call self%PureEFTwDE%init   ( 'EFTw'     , 'w'            )
        call self%PureEFTGamma1%init( 'EFTGamma1', '\gamma^{(1)}' )
        call self%PureEFTGamma2%init( 'EFTGamma2', '\gamma^{(2)}' )
        call self%PureEFTGamma3%init( 'EFTGamma3', '\gamma^{(3)}' )
        call self%PureEFTGamma4%init( 'EFTGamma4', '\gamma^{(4)}' )
        call self%PureEFTGamma5%init( 'EFTGamma5', '\gamma^{(5)}' )
        call self%PureEFTGamma6%init( 'EFTGamma6', '\gamma^{(6)}' )

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
        call self%PureEFTGamma4%init_from_file( Ini )
        call self%PureEFTGamma5%init_from_file( Ini )
        call self%PureEFTGamma6%init_from_file( Ini )

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
        self%parameter_number = self%parameter_number +self%PureEFTGamma4%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTGamma5%parameter_number
        self%parameter_number = self%parameter_number +self%PureEFTGamma6%parameter_number

    end subroutine EFTCAMBPureEFTstdComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBPureEFTstdFeedback( self )

        implicit none

        class(EFTCAMB_std_pure_EFT)  :: self   !< the base class

        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
        ! print model functions informations:
        write(*,*)
        if ( self%PureEFTmodelOmega  /= 0 ) write(*,'(a,I3)') '   PureEFTmodelOmega   =', self%PureEFTmodelOmega
        if ( self%EFTwDE             /= 0 ) write(*,'(a,I3)') '   EFTwDE              =', self%EFTwDE
        if ( self%PureEFTmodelGamma1 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma1  =', self%PureEFTmodelGamma1
        if ( self%PureEFTmodelGamma2 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma2  =', self%PureEFTmodelGamma2
        if ( self%PureEFTmodelGamma3 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma3  =', self%PureEFTmodelGamma3
        if ( self%PureEFTmodelGamma4 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma4  =', self%PureEFTmodelGamma4
        if ( self%PureEFTmodelGamma5 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma5  =', self%PureEFTmodelGamma5
        if ( self%PureEFTmodelGamma6 /= 0 ) write(*,'(a,I3)') '   PureEFTmodelGamma6  =', self%PureEFTmodelGamma6
        write(*,*)

        call self%PureEFTOmega%feedback()
        call self%PureEFTwDE%feedback()
        call self%PureEFTGamma1%feedback()
        call self%PureEFTGamma2%feedback()
        call self%PureEFTGamma3%feedback()
        call self%PureEFTGamma4%feedback()
        call self%PureEFTGamma5%feedback()
        call self%PureEFTGamma6%feedback()

    end subroutine EFTCAMBPureEFTstdFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

    end subroutine EFTCAMBPureEFTstdParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

    end subroutine EFTCAMBPureEFTstdParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

    end subroutine EFTCAMBPureEFTstdParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBPureEFTstdBackgroundEFTFunctions( self, a, &
        & EFTOmegaV, EFTOmegaP    , EFTOmegaPP, EFTOmegaPPP, &
        & EFTc     , EFTcdot      , &
        & EFTLambda, EFTLambdadot   &
        & )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self    !< the base class
        real(dl), intent(in)  :: a             !< the input scale factor
        real(dl), intent(out) :: EFTOmegaV     !< the value of Omega
        real(dl), intent(out) :: EFTOmegaP     !< the value of the derivative wrt scale factor of Omega
        real(dl), intent(out) :: EFTOmegaPP    !< the value of the second derivative wrt scale factor of Omega
        real(dl), intent(out) :: EFTOmegaPPP   !< the value of the third derivative wrt scale factor of Omega
        real(dl), intent(out) :: EFTc          !< the value of c*a^2/m_0^2
        real(dl), intent(out) :: EFTcdot       !< the value of \dot{c}*a^2/m_0^2. Derivative of c wrt conformal time
        real(dl), intent(out) :: EFTLambda     !< the value of \Lambda*a^2/m_0^2
        real(dl), intent(out) :: EFTLambdadot  !< the value of \dot{\Lambda}*a^2/m_0^2. Derivative of \Lambda wrt conformal time

        EFTOmegaV    = self%PureEFTOmega%value(a)
        EFTOmegaP    = self%PureEFTOmega%first_derivative(a)
        EFTOmegaPP   = self%PureEFTOmega%second_derivative(a)
        EFTOmegaPPP  = self%PureEFTOmega%third_derivative(a)
        EFTc         = 0._dl
        EFTcdot      = 0._dl
        EFTLambda    = 0._dl
        EFTLambdadot = 0._dl

    end subroutine EFTCAMBPureEFTstdBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBPureEFTstdSecondOrderEFTFunctions( self, a, &
        & EFTGamma1V, EFTGamma1P, &
        & EFTGamma2V, EFTGamma2P, &
        & EFTGamma3V, EFTGamma3P, &
        & EFTGamma4V, EFTGamma4P, EFTGamma4PP, &
        & EFTGamma5V, EFTGamma5P, &
        & EFTGamma6V, EFTGamma6P  &
        & )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self    !< the base class
        real(dl), intent(in)  :: a             !< the input scale factor
        real(dl), intent(out) :: EFTGamma1V    !< the value of Gamma 1
        real(dl), intent(out) :: EFTGamma1P    !< the value of the derivative wrt scale factor of Gamma 1
        real(dl), intent(out) :: EFTGamma2V    !< the value of Gamma 2
        real(dl), intent(out) :: EFTGamma2P    !< the value of the derivative wrt scale factor of Gamma 2
        real(dl), intent(out) :: EFTGamma3V    !< the value of Gamma 3
        real(dl), intent(out) :: EFTGamma3P    !< the value of the derivative wrt scale factor of Gamma 3
        real(dl), intent(out) :: EFTGamma4V    !< the value of Gamma 4
        real(dl), intent(out) :: EFTGamma4P    !< the value of the derivative wrt scale factor of Gamma 4
        real(dl), intent(out) :: EFTGamma4PP   !< the value of the second derivative wrt scale factor of Gamma 4
        real(dl), intent(out) :: EFTGamma5V    !< the value of Gamma 5
        real(dl), intent(out) :: EFTGamma5P    !< the value of the derivative wrt scale factor of Gamma 5
        real(dl), intent(out) :: EFTGamma6V    !< the value of Gamma 6
        real(dl), intent(out) :: EFTGamma6P    !< the value of the derivative wrt scale factor of Gamma 6

        EFTGamma1V  = self%PureEFTGamma1%value(a)
        EFTGamma1P  = self%PureEFTGamma1%first_derivative(a)
        EFTGamma2V  = self%PureEFTGamma2%value(a)
        EFTGamma2P  = self%PureEFTGamma2%first_derivative(a)
        EFTGamma3V  = self%PureEFTGamma3%value(a)
        EFTGamma3P  = self%PureEFTGamma3%first_derivative(a)
        EFTGamma4V  = self%PureEFTGamma4%value(a)
        EFTGamma4P  = self%PureEFTGamma4%first_derivative(a)
        EFTGamma4PP = self%PureEFTGamma4%second_derivative(a)
        EFTGamma5V  = self%PureEFTGamma5%value(a)
        EFTGamma5P  = self%PureEFTGamma5%first_derivative(a)
        EFTGamma6V  = self%PureEFTGamma6%value(a)
        EFTGamma6P  = self%PureEFTGamma6%first_derivative(a)

    end subroutine EFTCAMBPureEFTstdSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBPureEFTstdComputeDtauda( self, a, grhoa2, &
        & grhok, grhov, &
        & grhoc, grhob, &
        & grhog, grhornomass )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self                     !< the base class
        real(dl), intent(in)  :: a                              !< the input scale factor
        real(dl), intent(in)  :: grhoa2                         !< the input value of 8 \piG \rho_tot a^2
        real(dl), intent(in)  :: grhok                          !< the input value of curvature density
        real(dl), intent(in)  :: grhov                          !< the input value of DE density
        real(dl), intent(in)  :: grhoc                          !< the input value of CDM density
        real(dl), intent(in)  :: grhob                          !< the input value of Baryon density
        real(dl), intent(in)  :: grhog                          !< the input value of Radiation density
        real(dl), intent(in)  :: grhornomass                    !< the input value of massless neutrinos density

        real(dl)              :: EFTCAMBPureEFTstdComputeDtauda !< the output dtauda

        real(dl)              :: temp

        temp = grhoa2 +grhov*a*a*self%PureEFTwDE%integral(a)
        EFTCAMBPureEFTstdComputeDtauda = sqrt(3/temp)

    end function EFTCAMBPureEFTstdComputeDtauda

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_pure_EFT_std

!----------------------------------------------------------------------------------------

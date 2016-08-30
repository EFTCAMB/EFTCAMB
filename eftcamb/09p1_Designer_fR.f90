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

!> @file 09p1_Designer_fR.f90
!! This file contains the relevant code for designer f(R) models.


!----------------------------------------------------------------------------------------
!> This module contains the relevant code for designer f(R) models.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_designer_fR

    use precision
    use IniFile
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_neutral_parametrization_1D
    use EFTCAMB_constant_parametrization_1D
    use EFTCAMB_linear_parametrizations_1D
    use EFTCAMB_abstract_model_designer

    implicit none

    private

    public EFTCAMB_fR_designer

    !----------------------------------------------------------------------------------------
    !> This is the designer f(R) model.
    type, extends ( EFTCAMB_designer_model ) :: EFTCAMB_fR_designer

        ! theory parameters:
        real(dl) :: B0                                                    !< The present day value of B0.

        ! the pure EFT functions model selection flags:
        integer  :: EFTwDE                                                !< Model selection flag for designer f(R) w DE.

        ! the pure EFT functions:
        class( parametrized_function_1D ), allocatable :: PureEFTwDE      !< The pure EFT function w_DE.

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBDesignerFRReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBDesignerFRAllocateModelSelection      !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBDesignerFRInitModelParameters         !< subroutine taht initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBDesignerFRInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBDesignerFRInitBackground   !< subroutine that initializes the background of designer f(R).

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBDesignerFRComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBDesignerFRFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBDesignerFRParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBDesignerFRParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBDesignerFRParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBDesignerFRBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBDesignerFRSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_dtauda                    => EFTCAMBDesignerFRComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure :: compute_adotoa                    => EFTCAMBDesignerFRComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBDesignerFRComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.

    end type EFTCAMB_fR_designer

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBDesignerFRReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_fR_designer)  :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        ! read model selection flags:
        self%EFTwDE             = Ini_Read_Int_File( Ini, 'EFTwDE'             , 0 )

    end subroutine EFTCAMBDesignerFRReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBDesignerFRAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_fR_designer)  :: self   !< the base class

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

        ! initialize the names:
        call self%PureEFTwDE%set_name( 'EFTw'     , 'w'            )

    end subroutine EFTCAMBDesignerFRAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine taht initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBDesignerFRInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_fR_designer)  :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.

        stop 'IW'

    end subroutine EFTCAMBDesignerFRInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBDesignerFRInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_fR_designer)  :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        ! read w_DE parameters:
        call self%PureEFTwDE%init_from_file( Ini )
        ! read B0:
        self%B0 = Ini_Read_Double_File( Ini, 'EFTB0', 0._dl )

    end subroutine EFTCAMBDesignerFRInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of designer f(R).
    subroutine EFTCAMBDesignerFRInitBackground( self, params_cache )

        implicit none

        class(EFTCAMB_fR_designer)                   :: self          !< the base class
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters









        write(*,*) 'ciao Marco', params_cache%h0
        stop

    end subroutine EFTCAMBDesignerFRInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBDesignerFRComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_fR_designer)  :: self   !< the base class

        self%parameter_number = 1
        self%parameter_number = self%parameter_number +self%PureEFTwDE%parameter_number

    end subroutine EFTCAMBDesignerFRComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBDesignerFRFeedback( self )

        implicit none

        class(EFTCAMB_fR_designer)  :: self   !< the base class

        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
        ! print model functions informations:
        if ( self%EFTwDE /= 0 ) then
            write(*,*)
            write(*,'(a,I3)')  '   EFTwDE              =', self%EFTwDE
        end if
        write(*,*)
        write(*,'(a24,F12.6)') '   B0                  =', self%B0

        call self%PureEFTwDE%feedback()

    end subroutine EFTCAMBDesignerFRFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBDesignerFRParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_fR_designer) :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        stop 'IW: EFTCAMBDesignerFRParameterNames'

    end subroutine EFTCAMBDesignerFRParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBDesignerFRParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_fR_designer) :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

        stop 'IW: EFTCAMBDesignerFRParameterNamesLatex'

    end subroutine EFTCAMBDesignerFRParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBDesignerFRParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_fR_designer) :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

        stop 'IW: EFTCAMBDesignerFRParameterValues'

    end subroutine EFTCAMBDesignerFRParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBDesignerFRBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_designer)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%EFTOmegaV    = 0._dl
        eft_cache%EFTOmegaP    = 0._dl
        eft_cache%EFTOmegaPP   = 0._dl
        eft_cache%EFTOmegaPPP  = 0._dl
        eft_cache%EFTc         = 0._dl
        eft_cache%EFTLambda    = 0._dl
        eft_cache%EFTcdot      = 0._dl
        eft_cache%EFTLambdadot = 0._dl

    end subroutine EFTCAMBDesignerFRBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBDesignerFRSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_designer)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%EFTGamma1V  = 0._dl
        eft_cache%EFTGamma1P  = 0._dl
        eft_cache%EFTGamma2V  = 0._dl
        eft_cache%EFTGamma2P  = 0._dl
        eft_cache%EFTGamma3V  = 0._dl
        eft_cache%EFTGamma3P  = 0._dl
        eft_cache%EFTGamma4V  = 0._dl
        eft_cache%EFTGamma4P  = 0._dl
        eft_cache%EFTGamma4PP = 0._dl
        eft_cache%EFTGamma5V  = 0._dl
        eft_cache%EFTGamma5P  = 0._dl
        eft_cache%EFTGamma6V  = 0._dl
        eft_cache%EFTGamma6P  = 0._dl

    end subroutine EFTCAMBDesignerFRSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBDesignerFRComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_designer)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBDesignerFRComputeDtauda               !< the output dtauda

        real(dl) :: temp

        temp = eft_cache%grhoa2 +eft_par_cache%grhov*a*a*self%PureEFTwDE%integral(a)
        EFTCAMBDesignerFRComputeDtauda = sqrt(3/temp)

    end function EFTCAMBDesignerFRComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBDesignerFRComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_designer)                  :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%grhov_t = eft_par_cache%grhov*self%PureEFTwDE%integral(a)
        eft_cache%adotoa  = sqrt( ( eft_cache%grhom_t +eft_cache%grhov_t )/3._dl )

    end subroutine EFTCAMBDesignerFRComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBDesignerFRComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_designer)                  :: self          !< the base class
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

    end subroutine EFTCAMBDesignerFRComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_designer_fR

!----------------------------------------------------------------------------------------

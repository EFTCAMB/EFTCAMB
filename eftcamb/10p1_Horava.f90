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

!> @file 10p1_Horava.f90
!! This file contains the definition of low energy (LE) Horava gravity.
!! Please refer to the numerical notes for details.


!----------------------------------------------------------------------------------------
!> This module contains the definition of low energy (LE) Horava gravity.
!! Please refer to the numerical notes for details.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_LE_Horava

    use precision
    use IniFile
    use AMLutils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_model_designer

    implicit none

    private

    public EFTCAMB_Horava

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of low energy Horava gravity.
    type, extends ( EFTCAMB_designer_model ) :: EFTCAMB_Horava

        ! the model selection flag:
        logical   :: HoravaSolarSystem              !< Selects whether to consider the version of Horava gravity that evades solar system constraints.
                                                    !! Notice that when this is used the parameter Horava_xi is ignored.

        ! the model parameters:
        real(dl)  :: Horava_eta      !< first Horava model parameter.
        real(dl)  :: Horava_xi       !< second Horava model parameter.
        real(dl)  :: Horava_lambda   !< third Horava model parameter.

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBHoravaReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBHoravaAllocateModelSelection      !< subroutine that allocates the model selection. For Horava this is a dummy procedure.
        procedure :: init_model_parameters           => EFTCAMBHoravaInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBHoravaInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.
        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBHoravaComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBHoravaFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBHoravaParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBHoravaParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBHoravaParameterValues            !< subroutine that returns the i-th parameter value.
        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBHoravaBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBHoravaSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_dtauda                    => EFTCAMBHoravaComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure :: compute_adotoa                    => EFTCAMBHoravaComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBHoravaComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.
        procedure :: compute_pi_factors                => EFTCAMBHoravaComputePiFactors         !< subroutine that computes the pi field equations factors. For details refer to the numerical notes. Notice that in Horava we override it for reasons of numerical stability.
        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBHoravaAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_Horava

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBHoravaReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Horava)       :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        ! read model selection flags:
        self%HoravaSolarSystem = Ini_Read_Logical_File( Ini, 'HoravaSolarSystem', .False. )

    end subroutine EFTCAMBHoravaReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection. For Horava this is a dummy procedure.
    subroutine EFTCAMBHoravaAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_Horava) :: self !< the base class

    end subroutine EFTCAMBHoravaAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBHoravaInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_Horava)                                  :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

        real(dl), allocatable, dimension(:)                    :: temp
        integer :: num_params_function, num_params_temp, i

        if ( self%HoravaSolarSystem ) then
            self%Horava_eta    = array(1)
            self%Horava_lambda = array(2)
            self%Horava_xi     = 0.5_dl*self%Horava_eta
        else
            self%Horava_eta    = array(1)
            self%Horava_lambda = array(2)
            self%Horava_xi     = array(3)
        end if

    end subroutine EFTCAMBHoravaInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBHoravaInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Horava)          :: self   !< the base class
        type(TIniFile)                 :: Ini    !< Input ini file

        self%Horava_eta    = Ini_Read_Double_File( Ini, 'Horava_eta'   , 0._dl )
        self%Horava_lambda = Ini_Read_Double_File( Ini, 'Horava_lambda', 0._dl )
        if ( self%HoravaSolarSystem ) then
            self%Horava_xi = 0.5_dl*self%Horava_eta
        else
            self%Horava_xi = Ini_Read_Double_File( Ini, 'Horava_xi'    , 0._dl )
        end if

    end subroutine EFTCAMBHoravaInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBHoravaComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_Horava)  :: self   !< the base class

        if ( self%HoravaSolarSystem ) then
            self%parameter_number = 2
        else
            self%parameter_number = 3
        end if

    end subroutine EFTCAMBHoravaComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBHoravaFeedback( self )

        implicit none

        class(EFTCAMB_Horava)  :: self   !< the base class

        ! print general model informations:
        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number

        if ( self%HoravaSolarSystem ) then
            write(*,'(a)')  '   Horava with solar system constraints'
        end if

        ! print the values of the parameters:
        write(*,*)
        write(*,'(a23,a,F12.6)') '   Horava_eta          ', '=', self%Horava_eta
        write(*,'(a23,a,F12.6)') '   Horava_lambda       ', '=', self%Horava_lambda
        if ( .not. self%HoravaSolarSystem ) then
            write(*,'(a23,a,F12.6)') '   Horava_xi           ', '=', self%Horava_xi
        end if

    end subroutine EFTCAMBHoravaFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBHoravaParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_Horava)       :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            name = 'Horava_eta'
            return
        end if
        if ( i==2 ) then
            name = 'Horava_lambda'
            return
        end if
        if ( i==3 .and. self%HoravaSolarSystem ) then
            name = 'Horava_xi'
            return
        end if

    end subroutine EFTCAMBHoravaParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBHoravaParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_Horava)       :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names_latex.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            latexname = '\eta_{\rm Ho\v rava}'
            return
        end if
        if ( i==2 ) then
            latexname = '\lambda_{\rm Ho\v rava}'
            return
        end if
        if ( i==3 .and. self%HoravaSolarSystem ) then
            latexname = '\xi_{\rm Ho\v rava}'
            return
        end if

    end subroutine EFTCAMBHoravaParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBHoravaParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_Horava)       :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_value.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            value = self%Horava_eta
            return
        end if
        if ( i==2 ) then
            value = self%Horava_lambda
            return
        end if
        if ( i==3 .and. self%HoravaSolarSystem ) then
            value = self%Horava_xi
            return
        end if

    end subroutine EFTCAMBHoravaParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBHoravaBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Horava)                        :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        ! compute the EFT functions:
        eft_cache%EFTOmegaV    = self%Horava_eta/( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )
        eft_cache%EFTOmegaP    = 0._dl
        eft_cache%EFTOmegaPP   = 0._dl
        eft_cache%EFTOmegaPPP  = 0._dl
        eft_cache%EFTc         = -( 2._dl*self%Horava_xi -3._dl*self%Horava_lambda )*( eft_cache%Hdot -eft_cache%adotoa**2 )/( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )
        eft_cache%EFTLambda    = -eft_cache%grhov_t + 2._dl*( 3._dl*self%Horava_lambda - 2._dl*self%Horava_xi)*&
            &( 0.5_dl*eft_cache%adotoa**2 + eft_cache%Hdot)/( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )
        eft_cache%EFTcdot      = ( 3._dl*self%Horava_lambda - 2._dl*self%Horava_xi )*&
            &( eft_cache%Hdotdot -4._dl*eft_cache%adotoa*eft_cache%Hdot + 2._dl*eft_cache%adotoa**3 )/( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )
        eft_cache%EFTLambdadot = +2._dl*( 3._dl*self%Horava_lambda -2._dl*self%Horava_xi)*&
            &( eft_cache%Hdotdot -eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3 )/( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )

    end subroutine EFTCAMBHoravaBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBHoravaSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Horava)                        :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        ! compute the EFT functions:
        eft_cache%EFTGamma1V  = 0._dl
        eft_cache%EFTGamma1P  = 0._dl
        eft_cache%EFTGamma2V  = 0._dl
        eft_cache%EFTGamma2P  = 0._dl
        eft_cache%EFTGamma3V  = 2._dl*(self%Horava_lambda -self%Horava_xi)/( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )
        eft_cache%EFTGamma3P  = 0._dl
        eft_cache%EFTGamma4V  = 2._dl*self%Horava_xi/( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )
        eft_cache%EFTGamma4P  = 0._dl
        eft_cache%EFTGamma4PP = 0._dl
        eft_cache%EFTGamma5V  = 0._dl
        eft_cache%EFTGamma5P  = 0._dl
        eft_cache%EFTGamma6V  = self%Horava_eta/4._dl/( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )
        eft_cache%EFTGamma6P  = 0._dl

    end subroutine EFTCAMBHoravaSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBHoravaComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Horava)                        :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBHoravaComputeDtauda                        !< the output dtauda

        real(dl) :: temp

        temp = eft_cache%grhoa2 +a**4*( eft_par_cache%grhov &
            & +3._dl*eft_par_cache%h0_Mpc**2*( (self%Horava_eta +3._dl*self%Horava_lambda-2._dl*self%Horava_xi)/(2._dl*self%Horava_xi+2._dl-self%Horava_eta)))
        EFTCAMBHoravaComputeDtauda = sqrt(3/temp)/sqrt(( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )/(3._dl*self%Horava_lambda +2._dl))

    end function EFTCAMBHoravaComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBHoravaComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Horava)                        :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%grhov_t = ( eft_par_cache%grhov &
            & +3._dl*eft_par_cache%h0_Mpc**2*((self%Horava_eta +3._dl*self%Horava_lambda -2._dl*self%Horava_xi)/(2._dl*self%Horava_xi +2._dl -self%Horava_eta)))*a**2
        eft_cache%adotoa  = sqrt( ( eft_cache%grhom_t +eft_cache%grhov_t )/3._dl )
        eft_cache%adotoa  = eft_cache%adotoa*sqrt(( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )/(3._dl*self%Horava_lambda +2._dl))

    end subroutine EFTCAMBHoravaComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBHoravaComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Horava)                        :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%gpiv_t  = - eft_cache%grhov_t
        eft_cache%Hdot    = -0.5_dl*( eft_cache%adotoa**2 +eft_cache%gpresm_t +eft_cache%gpiv_t )

        ! IW POSSIBLE BUG
        !eft_cache%Hdotdot = eft_cache%adotoa*( ( eft_cache%grhob_t +eft_cache%grhoc_t)/6._dl +2._dl*( eft_cache%grhor_t +eft_cache%grhog_t)/3._dl ) &
        !    & +eft_cache%adotoa*eft_cache%grhov_t*( 1._dl/6._dl +self%Horava_wDE%value(a) +1.5_dl*self%Horava_wDE%value(a)**2 -0.5_dl*a*self%Horava_wDE%first_derivative(a) ) &
        !    & +eft_cache%adotoa*eft_cache%grhonu_tot/6._dl -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot

        eft_cache%Hdotdot = 2._dl*eft_cache%adotoa*eft_cache%Hdot &
            & + 0.5_dl*eft_cache%adotoa*(eft_cache%grhob_t + eft_cache%grhoc_t + 8._dl*(eft_cache%grhog_t+eft_cache%grhor_t)/3._dl)&
            & + eft_cache%adotoa/6._dl*eft_cache%grhonu_tot -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot

        eft_cache%Hdot     = eft_cache%Hdot*sqrt(( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )/(3._dl*self%Horava_lambda +2._dl))
        eft_cache%Hdotdot  = eft_cache%Hdotdot*sqrt(( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )/(3._dl*self%Horava_lambda +2._dl))

    end subroutine EFTCAMBHoravaComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the pi field equations factors. For details refer to the numerical notes. Notice that in Horava we override it for reasons of numerical stability.
    subroutine EFTCAMBHoravaComputePiFactors( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Horava)                        :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%EFTpiA1  = 0._dl
        eft_cache%EFTpiA2  = 0.5_dl*self%Horava_eta
        eft_cache%EFTpiB1  = 0._dl
        eft_cache%EFTpiB2  = eft_cache%adotoa*self%Horava_eta
        eft_cache%EFTpiC   = 0._dl
        eft_cache%EFTpiD1  = 0.5_dl*( 3._dl*self%Horava_lambda -2._dl*self%Horava_xi)*( eft_cache%adotoa**2 -eft_cache%Hdot ) + 0.5_dl*self%Horava_eta*( eft_cache%adotoa**2 +eft_cache%Hdot )
        eft_cache%EFTpiD2  = 0.5_dl*self%Horava_lambda*(1._dl+self%Horava_xi)
        eft_cache%EFTpiE   = 0.5_dl*eft_cache%k**2*self%Horava_lambda*( 1._dl +self%Horava_xi )*eft_cache%k*eft_cache%z + 0.5_dl*eft_cache%k*self%Horava_xi*eft_cache%dgq

    end subroutine EFTCAMBHoravaComputePiFactors

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBHoravaAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Horava)                        :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBHoravaAdditionalModelStability              !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBHoravaAdditionalModelStability = .True.

        if ( self%Horava_lambda > -2._dl/3._dl .and. self%Horava_lambda < 0._dl ) then
            EFTCAMBHoravaAdditionalModelStability = .False.
        end if
        if ( self%Horava_eta < 0._dl .or. self%Horava_eta > 2._dl*self%Horava_xi +2._dl ) then
            EFTCAMBHoravaAdditionalModelStability = .False.
        end if

    end function EFTCAMBHoravaAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_LE_Horava

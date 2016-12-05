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

!> @file 10p2_Cubic_Galileon.f90
!! This file contains the definition of the Cubic Galileon model.
!! Please refer to the numerical notes for details.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Cubic Galileon model.
!! Please refer to the numerical notes for details.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_full_Cubic_Galileon

    use precision
    use IniFile
    use AMLutils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_model_full
    ! use EFTCAMB_bilinear_parametrizations_2D

    implicit none

    private

    public EFTCAMB_Cubic_Galileon

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of the Cubic Galileon model.
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_Cubic_Galileon

        ! the model selection flag:
        logical   :: CubicGalileon         !< Selects Cubic Galileon model.

        ! the model parameters:
        real(dl)  :: CubicGalileon_c3      !< Cubic Galileon model parameter.

        ! ! the model functions:
        ! class( parametrized_function_2D ), allocatable :: G2    !< The Galileon function G2.
        ! class( parametrized_function_2D ), allocatable :: G3    !< The Galileon function G3.

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBCubicGalileonReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBCubicGalileonAllocateModelSelection      !< subroutine that allocates the model selection. For Horava this is a dummy procedure.
        procedure :: init_model_parameters           => EFTCAMBCubicGalileonInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBCubicGalileonInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.
        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBCubicGalileonComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBCubicGalileonFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBCubicGalileonParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBCubicGalileonParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBCubicGalileonParameterValues            !< subroutine that returns the i-th parameter value.
        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBCubicGalileonBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBCubicGalileonSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_dtauda                    => EFTCAMBCubicGalileonComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure :: compute_adotoa                    => EFTCAMBCubicGalileonComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBCubicGalileonComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.
        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBCubicGalileonAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_Cubic_Galileon

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBCubicGalileonReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Cubic_Galileon)       :: self   !< the base class
        type(TIniFile)                      :: Ini    !< Input ini file

        ! read model selection flags:
        self%CubicGalileon = Ini_Read_Logical_File( Ini, 'CubicGalileon', .False. )

    end subroutine EFTCAMBCubicGalileonReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBCubicGalileonAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_Cubic_Galileon) :: self !< the base class

        ! allocate G2 and G3:
        ! if ( allocated(self%G2) ) deallocate(self%G2)
        ! if ( allocated(self%G3) ) deallocate(self%G3)
        ! if (self%CubicGalileon) then
        !   allocate( bilinear_parametrization_2D:: self%G2 )
        !   call self%G2%set_param_names(['G2_phi  ', 'G2_X'], ['G2_{\phi}', 'G2_{\Chi}'])
        !   allocate( bilinear_parametrization_2D:: self%G3 )
        !   call self%G3%set_param_names(['G3_phi  ', 'G3_X'], ['G3_{\phi}', 'G3_{\Chi}'])
        ! else
        !   write(*,'(a,I3)') 'No model corresponding to CubicGalileon =', self%CubicGalileon
        !   write(*,'(a)')    'Please select an appropriate model.'
        ! end if

    end subroutine EFTCAMBCubicGalileonAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBCubicGalileonInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_Cubic_Galileon)                          :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

        real(dl), allocatable, dimension(:)                    :: temp

        self%CubicGalileon_c3    = array(1)


    end subroutine EFTCAMBCubicGalileonInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBCubicGalileonInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Cubic_Galileon)  :: self   !< the base class
        type(TIniFile)                 :: Ini    !< Input ini file

        self%CubicGalileon_c3    = Ini_Read_Double_File( Ini, 'CubicGalileon_c3'   , 0._dl )

    end subroutine EFTCAMBCubicGalileonInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBCubicGalileonComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_Cubic_Galileon)  :: self   !< the base class

        self%parameter_number = 1

    end subroutine EFTCAMBCubicGalileonComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBCubicGalileonFeedback( self, print_params )

        implicit none

        class(EFTCAMB_Cubic_Galileon)  :: self         !< the base class
        logical, optional              :: print_params !< optional flag that decised whether to print numerical values
                                                       !! of the parameters.

        logical                        :: print_params_temp

        ! print general model informations:
        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number

        ! print the values of the parameters:
        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

        if ( print_params_temp ) then
            write(*,*)
            write(*,'(a23,a,F12.6)') '   CubicGalileon_c3    ', '=', self%CubicGalileon_c3
        end if

    end subroutine EFTCAMBCubicGalileonFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBCubicGalileonParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_Cubic_Galileon) :: self   !< the base class
        integer     , intent(in)      :: i      !< the index of the parameter
        character(*), intent(out)     :: name   !< the output name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            name = 'CubicGalileon_c3'
            return
        end if

    end subroutine EFTCAMBCubicGalileonParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBCubicGalileonParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_Cubic_Galileon) :: self       !< the base class
        integer     , intent(in)      :: i          !< The index of the parameter
        character(*), intent(out)     :: latexname  !< the output latex name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names_latex.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            latexname = 'c_3'
            return
        end if

    end subroutine EFTCAMBCubicGalileonParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBCubicGalileonParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_Cubic_Galileon) :: self   !< the base class
        integer , intent(in)          :: i      !< The index of the parameter
        real(dl), intent(out)         :: value  !< the output value of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_value.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            value = self%CubicGalileon_c3
            return
        end if

    end subroutine EFTCAMBCubicGalileonParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBCubicGalileonBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Cubic_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        ! compute the EFT functions:
        eft_cache%EFTOmegaV    = 0._dl
        eft_cache%EFTOmegaP    = 0._dl
        eft_cache%EFTOmegaPP   = 0._dl
        eft_cache%EFTOmegaPPP  = 0._dl
        eft_cache%EFTc         = -0.5_dl*eft_cache%adotoa**2*eft_cache%phiP**2-0.5_dl*self%CubicGalileon_c3*eft_cache%adotoa**2*eft_cache%phiP**2*&
            &((3._dl*eft_cache%adotoa**2-eft_cache%Hdot)*eft_cache%phiP/a - eft_cache%adotoa**2*eft_cache%phiPP)
        eft_cache%EFTLambda    = -0.5_dl*eft_cache%adotoa**2*eft_cache%phiP**2+2._dl*self%CubicGalileon_c3*eft_cache%adotoa**2*eft_cache%phiP**2*&
            &( eft_cache%Hdot/a*eft_cache%phiP + eft_cache%Hdot**2*eft_cache%phiPP )
        !
        ! eft_cache%EFTcdot      =
        ! eft_cache%EFTLambdadot =

    end subroutine EFTCAMBCubicGalileonBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBCubicGalileonSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Cubic_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        !IW
        ! compute the EFT functions:
        eft_cache%EFTGamma1V  = 0.25_dl*self%CubicGalileon_c3*eft_cache%adotoa**2*eft_cache%phiP**2*(( 3._dl*eft_cache%adotoa**2 -eft_cache%Hdot)*eft_cache%phiP/a -&
            & eft_cache%adotoa**2*eft_cache%phiPP )
        eft_cache%EFTGamma1P  = 0._dl
        eft_cache%EFTGamma2V  = self%CubicGalileon_c3*eft_cache%adotoa**3*eft_cache%phiP**3
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

    end subroutine EFTCAMBCubicGalileonSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBCubicGalileonComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Cubic_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBCubicGalileonComputeDtauda                        !< the output dtauda

        real(dl) :: temp
        !IW
    end function EFTCAMBCubicGalileonComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBCubicGalileonComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Cubic_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        !IW
        ! eft_cache%grhov_t = ( eft_par_cache%grhov &
        !     & +3._dl*eft_par_cache%h0_Mpc**2*((self%Horava_eta +3._dl*self%Horava_lambda -2._dl*self%Horava_xi)/(2._dl*self%Horava_xi +2._dl -self%Horava_eta)))*a**2
        ! eft_cache%adotoa  = sqrt( ( eft_cache%grhom_t +eft_cache%grhov_t )/3._dl )
        ! eft_cache%adotoa  = eft_cache%adotoa*sqrt(( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )/(3._dl*self%Horava_lambda +2._dl))

    end subroutine EFTCAMBCubicGalileonComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBCubicGalileonComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Cubic_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        !IW
        ! eft_cache%gpiv_t   = -eft_cache%grhov_t
        ! eft_cache%Hdot     = -0.5_dl*( eft_cache%adotoa**2 +eft_cache%gpresm_t +eft_cache%gpiv_t )
        ! eft_cache%Hdotdot  = eft_cache%adotoa*( ( eft_cache%grhob_t +eft_cache%grhoc_t)/6._dl +2._dl*( eft_cache%grhor_t +eft_cache%grhog_t)/3._dl ) &
        !     & +2._dl*eft_cache%adotoa*eft_cache%grhov_t/3._dl &
        !     & +eft_cache%adotoa*eft_cache%grhonu_tot/6._dl -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot
        !
        ! eft_cache%Hdot     = eft_cache%Hdot*sqrt(( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )/(3._dl*self%Horava_lambda +2._dl))
        ! eft_cache%Hdotdot  = eft_cache%Hdotdot*sqrt(( 2._dl +2._dl*self%Horava_xi -self%Horava_eta )/(3._dl*self%Horava_lambda +2._dl))

    end subroutine EFTCAMBCubicGalileonComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBCubicGalileonAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Cubic_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBCubicGalileonAdditionalModelStability       !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBCubicGalileonAdditionalModelStability = .True.

        ! IW

    end function EFTCAMBCubicGalileonAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_full_Cubic_Galileon

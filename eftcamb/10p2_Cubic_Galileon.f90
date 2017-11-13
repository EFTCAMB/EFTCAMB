!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2017 by the EFTCAMB authors
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

!> @author Simone Peirone, Bin Hu, Marco Raveri

module EFTCAMB_full_Cubic_Galileon

    use precision
    use IniFile
    use AMLutils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_model_full

    implicit none

    private

    public EFTCAMB_Cubic_Galileon

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of the Cubic Galileon model.
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_Cubic_Galileon

        ! the model parameters:
        real(dl)  :: CubicGalileon_c2      !< Cubic Galileon model parameter \f$c_2\f$
        real(dl)  :: CubicGalileon_c3      !< Cubic Galileon model parameter \f$c_3\f$
        real(dl)  :: csi                   !< Cubic Galileon background parameter \f$\xi\f$ deriving from the tracker solution

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBCubicGalileonReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBCubicGalileonAllocateModelSelection      !< subroutine that allocates the model selection. For Horava this is a dummy procedure.
        procedure :: init_model_parameters           => EFTCAMBCubicGalileonInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBCubicGalileonInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBCubicGalileonInitBackground              !< subroutine that initializes the background of Cubic Galileon.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBCubicGalileonComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBCubicGalileonFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBCubicGalileonParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBCubicGalileonParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBCubicGalileonParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBCubicGalileonBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBCubicGalileonSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_adotoa                    => EFTCAMBCubicGalileonComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBCubicGalileonComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.

        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBCubicGalileonAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_Cubic_Galileon

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBCubicGalileonReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Cubic_Galileon)       :: self   !< the base class
        type(TIniFile)                      :: Ini    !< Input ini file

    end subroutine EFTCAMBCubicGalileonReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBCubicGalileonAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_Cubic_Galileon) :: self !< the base class

    end subroutine EFTCAMBCubicGalileonAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    !! Nothing needs to be done but procedure present because it is deferred.
    subroutine EFTCAMBCubicGalileonInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_Cubic_Galileon)                          :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

    end subroutine EFTCAMBCubicGalileonInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBCubicGalileonInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Cubic_Galileon)  :: self   !< the base class
        type(TIniFile)                 :: Ini    !< Input ini file

    end subroutine EFTCAMBCubicGalileonInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of Cubic Galileon.
    subroutine EFTCAMBCubicGalileonInitBackground( self, params_cache, feedback_level, success )

        implicit none

        class(EFTCAMB_Cubic_Galileon)                :: self           !< the base class
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                      , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                      , intent(out)   :: success        !< wether the background initialization succeded or not

        real(dl) :: Omega_phi0

        Omega_phi0 = params_cache%omegav

        ! SP: Cubic -> just c_3
        self%csi = sqrt( 6._dl*Omega_phi0 )
        self%CubicGalileon_c2 = -1.0_dl
        self%CubicGalileon_c3 =  1.0_dl/(6.0_dl*self%csi )
        call self%feedback()

        success=.true.

    end subroutine EFTCAMBCubicGalileonInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBCubicGalileonComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_Cubic_Galileon)  :: self   !< the base class

        self%parameter_number = 0

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
        if (self%CubicGalileon_c2 == -1._dl) then

          write(*,*)
          write(*,'(a,a)')    '   Model              =  ', self%name
          write(*,'(a,I3)')   '   Number of params   ='  , self%parameter_number
          write(*,'(a,F12.6)')   '                xi    ='  , self%csi
          write(*,'(a,F12.6)')   '                c2    ='  , self%CubicGalileon_c2
          write(*,'(a,F12.6)')   '                c3    ='  , self%CubicGalileon_c3

        end if

        ! print the values of the parameters:
        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
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
        if ( i==0 ) then
            name = 'no_name'
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
        if ( i==0 ) then
            latexname = 'noname'
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
        if ( i==0 ) then
            value = 0._dl
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

        real(dl)    :: psi, psiprime, psiprimeprime, a2, phip1, phip2, phip3, c3

        a2 = a*a

        if(a==0._dl) then
            return
        else if (eft_cache%adotoa==0._dl) then
            call self%compute_adotoa( a, eft_par_cache, eft_cache )
            call self%compute_H_derivs( a, eft_par_cache, eft_cache )
        end if

        ! compute the psi field and its derivatives
        psi           = self%csi*(eft_par_cache%h0_Mpc)**2*a2/eft_cache%adotoa**2
        psiprime      = 2._dl*self%csi*(eft_par_cache%h0_Mpc)**2*a/eft_cache%adotoa**4*( eft_cache%adotoa**2 -eft_cache%Hdot )
        psiprimeprime = 2._dl*self%csi*(eft_par_cache%h0_Mpc)**2/eft_cache%adotoa**4*( eft_cache%adotoa**2 -3._dl*eft_cache%Hdot &
            & +4._dl*( eft_cache%Hdot/eft_cache%adotoa )**2 -eft_cache%Hdotdot/eft_cache%adotoa )

        ! scalar field conversion
        phip1 = psi/a
        phip2 = ( psiprime*a -psi )/a2
        phip3 = ( psiprimeprime*a2 -2._dl*psiprime*a +2._dl*psi )/a**3

        c3= -1.0_dl*self%CubicGalileon_c3/(eft_par_cache%h0_Mpc)**2

        ! compute the background EFT functions:
        eft_cache%EFTOmegaV    = 0._dl
        eft_cache%EFTOmegaP    = 0._dl
        eft_cache%EFTOmegaPP   = 0._dl
        eft_cache%EFTOmegaPPP  = 0._dl
        eft_cache%EFTc         = 0.5_dl*self%CubicGalileon_c2*a2*eft_cache%adotoa**2*phip1**2 -c3*a2*eft_cache%adotoa**2*phip1**2*&
            & ( ( 3._dl*eft_cache%adotoa**2 -eft_cache%Hdot )*phip1/a -eft_cache%adotoa**2*phip2 )
        eft_cache%EFTLambda    = 0.5_dl*self%CubicGalileon_c2*a2*eft_cache%adotoa**2*phip1**2 +2._dl*c3*a2*eft_cache%adotoa**2*phip1**2*&
            & ( eft_cache%Hdot/a*phip1 +eft_cache%adotoa**2*phip2 )
        eft_cache%EFTcdot      = self%CubicGalileon_c2*a2*( eft_cache%adotoa*eft_cache%Hdot*phip1**2 +a*eft_cache%adotoa**3*phip1*phip2 ) &
            & -c3*a2*( 2._dl*( eft_cache%adotoa*eft_cache%Hdot*phip1**2 +a*eft_cache%adotoa**3*phip1*phip2 )*( ( 3._dl*eft_cache%adotoa**2 &
            & -eft_cache%Hdot )*phip1/a -eft_cache%adotoa**2*phip2 ) +eft_cache%adotoa**2*phip1**2*( ( 6._dl*eft_cache%adotoa*eft_cache%Hdot -eft_cache%Hdotdot )&
            & *phip1/a +eft_cache%adotoa*( 3._dl*eft_cache%adotoa**2 -eft_cache%Hdot )*( phip2 -phip1/a ) -2._dl*eft_cache%adotoa*eft_cache%Hdot*phip2 &
            & -a*eft_cache%adotoa**3*phip3 ) )
        eft_cache%EFTLambdadot = self%CubicGalileon_c2*a2*( eft_cache%adotoa*eft_cache%Hdot*phip1**2 +a*eft_cache%adotoa**3*phip1*phip2 ) &
            & +4._dl*c3*a2*( eft_cache%adotoa*eft_cache%Hdot*phip1**2 +a*eft_cache%adotoa**3*phip1*phip2 )*( eft_cache%Hdot*phip1/a &
            & +eft_cache%adotoa**2*phip2 ) +2._dl*c3*a2*eft_cache%adotoa**2*phip1**2*( eft_cache%Hdotdot*phip1/a +eft_cache%adotoa*eft_cache%Hdot*( 3._dl*phip2 &
            & -phip1/a ) +a*eft_cache%adotoa**3*phip3 )

    end subroutine EFTCAMBCubicGalileonBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the second order EFT functions at a given time.
    subroutine EFTCAMBCubicGalileonSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Cubic_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: psi, psiprime, psiprimeprime, a2, phip1, phip2, phip3, c3

        a2 = a*a

        if(a*eft_cache%adotoa==0._dl) return

        ! compute the psi field and its derivatives
        psi           = self%csi*(eft_par_cache%h0_Mpc)**2*a2/eft_cache%adotoa**2
        psiprime      = 2._dl*self%csi*(eft_par_cache%h0_Mpc)**2*a/eft_cache%adotoa**4*( eft_cache%adotoa**2 -eft_cache%Hdot )
        psiprimeprime = 2._dl*self%csi*(eft_par_cache%h0_Mpc)**2/eft_cache%adotoa**4*( eft_cache%adotoa**2 -3._dl*eft_cache%Hdot &
            & +4._dl*( eft_cache%Hdot/eft_cache%adotoa )**2 -eft_cache%Hdotdot/eft_cache%adotoa )

        ! scalar field conversion
        phip1 = psi/a
        phip2 = ( psiprime*a -psi )/a2
        phip3 = ( psiprimeprime*a2 -2._dl*psiprime*a +2._dl*psi )/a**3

        c3    = -1.0_dl*self%CubicGalileon_c3/(eft_par_cache%h0_Mpc)**2

        ! compute the second order EFT functions:
        eft_cache%EFTGamma1V  = -0.5_dl*c3/(eft_par_cache%h0_Mpc)**2*eft_cache%adotoa**2*phip1**2*( ( 3._dl*eft_cache%adotoa**2 &
            & +eft_cache%Hdot )*phip1/a +eft_cache%adotoa**2*phip2 )
        eft_cache%EFTGamma1P  = -0.5_dl*c3/(eft_par_cache%h0_Mpc)**2*( 2._dl*( eft_cache%Hdot*phip1**2/a +eft_cache%adotoa**2*phip1*phip2 )&
            & *( ( 3._dl*eft_cache%adotoa**2 +eft_cache%Hdot )*phip1/a +eft_cache%adotoa**2*phip2 ) + eft_cache%adotoa**2*phip1**2*( ( 6._dl*eft_cache%Hdot &
            & +eft_cache%Hdotdot/eft_cache%adotoa )*phip1/a2 +( 3._dl*eft_cache%adotoa**2 +eft_cache%Hdot )*( phip2/a -phip1/a2 ) +2._dl*eft_cache%Hdot/a*phip2&
            & +eft_cache%adotoa**2*phip3 ) )
        eft_cache%EFTGamma2V  = 2._dl*c3/(eft_par_cache%h0_Mpc)*( eft_cache%adotoa*phip1 )**3
        eft_cache%EFTGamma2P  = 6._dl*c3/(eft_par_cache%h0_Mpc)*eft_cache%adotoa*phip1**2*( eft_cache%Hdot/a*phip1 +eft_cache%adotoa**2*phip2 )
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
    !> Subroutine that computes adotoa = H.
    subroutine EFTCAMBCubicGalileonComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Cubic_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp, a2, Omega_tot
        integer     :: nu_i

        a2 = a*a

        Omega_tot = ( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-3) + ( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-4) +eft_cache%grhonu_tot/(3._dl*eft_par_cache%h0_Mpc**2*a2)

        temp = 0.5_dl*a2*(eft_par_cache%h0_Mpc)**2*( Omega_tot + sqrt( Omega_tot**2 +4._dl*eft_par_cache%omegav ) )
        eft_cache%adotoa = sqrt( temp )

    end subroutine EFTCAMBCubicGalileonComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBCubicGalileonComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Cubic_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp, a2, Omega_tot, Omega_tot_prime, Omega_tot_primeprime, Omega_phi0
        integer     :: nu_i

        a2 = a*a

        Omega_tot = ( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-3) + ( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-4) +eft_cache%grhonu_tot/(3._dl*eft_par_cache%h0_Mpc**2*a2)
        Omega_tot_prime = -3._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-4) -4._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-5) &
                          & -(eft_cache%grhonu_tot+eft_cache%gpinu_tot)/(eft_par_cache%h0_Mpc**2*a2*a)
        Omega_tot_primeprime = 12._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-5) +20._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-6)&
                          & +(4._dl*(eft_cache%grhonu_tot+eft_cache%gpinu_tot)-eft_cache%gpinudot_tot/eft_cache%adotoa )/(eft_par_cache%h0_Mpc**2*a2**2)

        Omega_phi0 = eft_par_cache%omegav


        eft_cache%Hdot = eft_cache%adotoa**2 +0.25_dl*(eft_par_cache%h0_Mpc)**2*a**3*( 1._dl + Omega_tot/sqrt( Omega_tot**2 +4._dl*Omega_phi0 ) )*Omega_tot_prime
        eft_cache%Hdotdot = 2._dl*eft_cache%adotoa*eft_cache%Hdot +3._dl*eft_cache%adotoa*( eft_cache%Hdot -eft_cache%adotoa**2 ) +0.25_dl*(eft_par_cache%h0_Mpc)**2*eft_cache%adotoa*a2**2&
            & *( ( 1._dl +Omega_tot/sqrt( Omega_tot**2 +4._dl*Omega_phi0 ) )*Omega_tot_primeprime +Omega_tot_prime**2&
            & *( 4._dl*Omega_phi0/( Omega_tot**2 +4._dl*Omega_phi0 )**( 1.5_dl ) ) )

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

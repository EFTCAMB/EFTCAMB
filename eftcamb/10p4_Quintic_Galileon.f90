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

!> @file 10p4_Quintic_Galileon.f90
!! This file contains the definition of the Quintic Galileon model.
!! Please refer to the numerical notes for details.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Quintic Galileon model.
!! Please refer to the numerical notes for details.

!> @author Simone Peirone, Bin Hu, Marco Raveri

module EFTCAMB_full_Quintic_Galileon

    use precision
    use IniFile
    use AMLutils
    use EFTCAMB_cache
    use EFT_def
    use equispaced_linear_interpolation_1D
    use EFTCAMB_abstract_model_full

    implicit none

    private

    public EFTCAMB_Quintic_Galileon

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of the Quintic Galileon model.
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_Quintic_Galileon

        ! the model parameters:
        real(dl)  :: QuinticGalileon_c2      !< Quintic Galileon model parameter \f$c_2\f$
        real(dl)  :: QuinticGalileon_c3      !< Quintic Galileon model parameter \f$c_3\f$
        real(dl)  :: QuinticGalileon_c4      !< Quintic Galileon model parameter \f$c_4\f$
        real(dl)  :: QuinticGalileon_c5      !< Quintic Galileon model parameter \f$c_5\f$
        real(dl)  :: csi                     !< Quintic Galileon background parameter \f$\xi\f$ deriving from the tracker solution

        ! the interpolated EFT functions that come out of the background sover:
        type(equispaced_linear_interpolate_function_1D) :: EFTOmega       !< The interpolated function Omega (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTLambda      !< The interpolated function Lambda (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTc           !< The interpolated function c (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTgamma1      !< The interpolated function \f$\gamma_1\f$ (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTgamma2      !< The interpolated function \f$\gamma_2\f$ (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTgamma3      !< The interpolated function \f$\gamma_3\f$ (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTgamma4      !< The interpolated function \f$\gamma_4\f$ (and derivatives).

        ! some designer parameters:
        integer  :: designer_num_points = 1000                             !< Number of points sampled by the designer code.
        real(dl) :: x_initial           = log(10._dl**(-10._dl))           !< log(a start)
        real(dl) :: x_final             = log(2._dl)                       !< log(a final)

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBQuinticGalileonReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBQuinticGalileonAllocateModelSelection      !< subroutine that allocates the model selection. For Horava this is a dummy procedure.
        procedure :: init_model_parameters           => EFTCAMBQuinticGalileonInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBQuinticGalileonInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBQuinticGalileonInitBackground              !< subroutine that initializes the background of Quintic Galileon.
        procedure :: solve_background                => EFTCAMBQuinticGalileonSolveBackground             !< subroutine that solves the background equations.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBQuinticGalileonComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBQuinticGalileonFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBQuinticGalileonParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBQuinticGalileonParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBQuinticGalileonParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBQuinticGalileonBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBQuinticGalileonSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_adotoa                    => EFTCAMBQuinticGalileonComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBQuinticGalileonComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.

        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBQuinticGalileonAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_Quintic_Galileon

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBQuinticGalileonReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Quintic_Galileon)         :: self   !< the base class
        type(TIniFile)                          :: Ini    !< Input ini file

    end subroutine EFTCAMBQuinticGalileonReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBQuinticGalileonAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_Quintic_Galileon) :: self !< the base class

    end subroutine EFTCAMBQuinticGalileonAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    !! Nothing needs to be done but procedure present because it is deferred.
    subroutine EFTCAMBQuinticGalileonInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_Quintic_Galileon)                            :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in)     :: array  !< input array with the values of the parameters of the model.

        self%csi                = array(1)
        self%QuinticGalileon_c3 = array(2)

    end subroutine EFTCAMBQuinticGalileonInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBQuinticGalileonInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Quintic_Galileon)  :: self   !< the base class
        type(TIniFile)                   :: Ini    !< Input ini file

        self%csi                = Ini_Read_Double_File( Ini, 'Quintic_Galileon_xi', 0._dl )
        self%QuinticGalileon_c3 = Ini_Read_Double_File( Ini, 'Quintic_Galileon_c3', 0._dl )

    end subroutine EFTCAMBQuinticGalileonInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of Quintic Galileon.
    subroutine EFTCAMBQuinticGalileonInitBackground( self, params_cache, feedback_level, success )

        implicit none

        class(EFTCAMB_Quintic_Galileon)              :: self           !< the base class
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                      , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                      , intent(out)   :: success        !< wether the background initialization succeded or not

        real(dl) :: Omega_phi0

        Omega_phi0 = params_cache%omegav

        ! Quintic -> all c_i
        self%QuinticGalileon_c2 = -1.0_dl
        self%QuinticGalileon_c5= 4._dl/3._dl*Omega_phi0*self%csi**(-5)+ 1._dl/3._dl*self%QuinticGalileon_c2*self%csi**(-3)+2._dl/3._dl*self%QuinticGalileon_c3*self%csi**(-2)
        self%QuinticGalileon_c4= -10._dl/9._dl*Omega_phi0*self%csi**(-4)- 1._dl/3._dl*self%QuinticGalileon_c2*self%csi**(-2)-8._dl/9._dl*self%QuinticGalileon_c3*self%csi**(-1)

        call self%feedback()

        ! initialize interpolating functions:
        call self%EFTOmega%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTLambda%initialize ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTc%initialize      ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTgamma1%initialize ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTgamma2%initialize ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTgamma3%initialize ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTgamma4%initialize ( self%designer_num_points, self%x_initial, self%x_final )

        ! solve the background equations and store the solution:
        call self%solve_background( params_cache, success=success )

        success = .true.

    end subroutine EFTCAMBQuinticGalileonInitBackground

    subroutine EFTCAMBQuinticGalileonSolveBackground( self, params_cache, success )

        implicit none

        class(EFTCAMB_Quintic_Galileon)              :: self          !< the base class.
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
        logical , intent(out)                        :: success       !< whether the calculation ended correctly or not
        real(dl) :: PPlus, yPlus, CoeffA_Part, yStar, x
        integer  :: i
        real(dl) :: t1, t2

        t1  = self%EFTOmega%x(1)
        call output(params_cache,  1, t1 )
        do i=1, self%EFTOmega%num_points-1

            t2 = self%EFTOmega%x(i+1)
            call output(params_cache,  i+1, t2)

        end do

        success =.true.
        return

    contains

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that takes the solution of the background and stores the values of the EFT functions.
      subroutine output( eft_par_cache, ind, x)

         implicit none

         type(EFTCAMB_parameter_cache), intent(in):: eft_par_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
         integer , intent(in)                     :: ind    !< index of the EFT functions interpolation tables to fill.
         real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed.

         real(dl) :: a,Omega_tot,Omega_tot_prime,Omega_tot_primeprime, Omega_tot_primeprimeprime,Omega_tot_primeprimeprimeprime
         real(dl) :: rhonu_tot, presnu_tot, presnudot_tot, presnudotdot_tot,presnudotdotdot_tot,Omega_phi0
         real(dl) :: rhonu, presnu, grhormass_t,presnudotdotdot,presnudotdot,presnudot, psi,psiprime, psiprimeprime, psiprimeprimeprime
         real(dl) :: adotoa, Hdot,Hdotdot,Hddd, Hdddd, psiprimeprimeprimeprime
         real(dl) :: phip1, phip2, phip3, phip4, phip5, m0, a2, c3
         integer  :: nu_i

         a = Exp(x)
         a2 = a*a
         m0=1._dl

         rhonu_tot  = 0._dl
         presnu_tot = 0._dl
         if ( eft_par_cache%Num_Nu_Massive /= 0 ) then
           do nu_i = 1, eft_par_cache%Nu_mass_eigenstates

             rhonu           = 0._dl
             presnu          = 0._dl
             grhormass_t = eft_par_cache%grhormass(nu_i)/a**2

             call eft_par_cache%Nu_background(a*eft_par_cache%nu_masses(nu_i),rhonu,presnu)
             rhonu_tot  = rhonu_tot + grhormass_t*rhonu
             presnu_tot = presnu_tot + grhormass_t*presnu

            end do
         end if

         Omega_tot = ( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-3) + ( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-4) +rhonu_tot/(3._dl*eft_par_cache%h0_Mpc**2*a2)
         adotoa = sqrt( 0.5_dl*a2*(eft_par_cache%h0_Mpc)**2*( Omega_tot + sqrt( Omega_tot**2 +4._dl*eft_par_cache%omegav ) ) )
         Omega_phi0 = eft_par_cache%omegav
         Omega_tot_prime = -3._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-4) -4._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-5) &
                                  & -(rhonu_tot+presnu_tot)/(eft_par_cache%h0_Mpc**2*a2*a)
         Hdot = adotoa**2 +0.25_dl*(eft_par_cache%h0_Mpc)**2*a**3*( 1._dl + Omega_tot/sqrt( Omega_tot**2 +4._dl*Omega_phi0 ) )*Omega_tot_prime

         rhonu_tot  = 0._dl
         presnu_tot = 0._dl
         presnudot_tot = 0._dl
         presnudotdot_tot = 0._dl
         if ( eft_par_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, eft_par_cache%Nu_mass_eigenstates

              rhonu           = 0._dl
              presnu          = 0._dl
              presnudot       = 0._dl
              presnudotdot    = 0._dl
              grhormass_t = eft_par_cache%grhormass(nu_i)/a**2

              call eft_par_cache%Nu_background(a*eft_par_cache%nu_masses(nu_i),rhonu,presnu)
              rhonu_tot  = rhonu_tot + grhormass_t*rhonu
              presnu_tot = presnu_tot + grhormass_t*presnu
              presnudot = params_cache%Nu_pidot(a*params_cache%nu_masses(nu_i),adotoa,presnu)
              presnudotdot = eft_par_cache%Nu_pidotdot(a*eft_par_cache%nu_masses(nu_i),adotoa,Hdot,presnu,presnudot)
              presnudot_tot  = presnudot_tot + grhormass_t*(presnudot -4._dl*adotoa*presnu)
              presnudotdot_tot = presnudotdot_tot + grhormass_t*(presnudotdot -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))

            end do
         end if

         Omega_tot_primeprime = 12._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-5) +20._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-6)&
                                  & +(4._dl*(rhonu_tot+presnu_tot)-presnudot_tot/adotoa )/(eft_par_cache%h0_Mpc**2*a2**2)

         Hdotdot = 2._dl*adotoa*Hdot +3._dl*adotoa*( Hdot -adotoa**2 ) +0.25_dl*(eft_par_cache%h0_Mpc)**2*adotoa*a2**2&
                       & *( ( 1._dl +Omega_tot/sqrt( Omega_tot**2 +4._dl*Omega_phi0 ) )*Omega_tot_primeprime +Omega_tot_prime**2&
                       & *( 4._dl*Omega_phi0/( Omega_tot**2 +4._dl*Omega_phi0 )**( 1.5_dl ) ) )
         !
         if ( a == 0._dl ) then
             return
         else if ( adotoa  == 0._dl ) then
             if  ( adotoa  == 0._dl ) return
             if  ( Hdot    == 0._dl ) return
             if  ( Hdotdot == 0._dl ) return
         end if
         rhonu_tot  = 0._dl
         presnu_tot = 0._dl
         presnudot_tot = 0._dl
         presnudotdot_tot = 0._dl
         presnudotdotdot_tot = 0._dl
         if ( eft_par_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, eft_par_cache%Nu_mass_eigenstates

              rhonu           = 0._dl
              presnu          = 0._dl
              presnudot       = 0._dl
              presnudotdot    = 0._dl
              presnudotdotdot = 0._dl
              grhormass_t = eft_par_cache%grhormass(nu_i)/a**2

              call eft_par_cache%Nu_background(a*eft_par_cache%nu_masses(nu_i),rhonu,presnu)
              presnudot = params_cache%Nu_pidot(a*params_cache%nu_masses(nu_i),adotoa,presnu)
              presnudotdot = eft_par_cache%Nu_pidotdot(a*eft_par_cache%nu_masses(nu_i),adotoa,Hdot,presnu,presnudot)
              presnudotdotdot = eft_par_cache%Nu_pidotdotdot(a*eft_par_cache%nu_masses(nu_i),adotoa,Hdot,Hdotdot,presnu,presnudot,presnudotdot )
              rhonu_tot  = rhonu_tot + grhormass_t*rhonu
              presnu_tot = presnu_tot + grhormass_t*presnu
              presnudot_tot  = presnudot_tot + grhormass_t*(presnudot -4._dl*adotoa*presnu)
              presnudotdot_tot = presnudotdot_tot + grhormass_t*(presnudotdot -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))
              presnudotdotdot_tot = presnudotdotdot_tot + grhormass_t*( presnudotdotdot -12._dl*adotoa*presnudotdot &
                  & -64._dl*adotoa**3*presnu -12._dl*Hdot*presnudot +48._dl*adotoa**2*presnudot -4._dl*Hdotdot*presnu +48._dl*adotoa*Hdot*presnu)

            end do
        end if

        Omega_tot_primeprimeprime = -60._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-6) -120._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-7)&
                    & +(-20._dl*(rhonu_tot+presnu_tot) + (6._dl/adotoa +Hdot/adotoa**3)*presnudot_tot  -1._dl/adotoa**2*presnudotdot_tot)/(eft_par_cache%h0_Mpc**2*a**5)

        Hddd = 9._dl*adotoa*Hdotdot -26._dl*adotoa**2*Hdot +Hdot*Hdotdot/adotoa &
                    &+12._dl*adotoa**4 +0.25_dl*(eft_par_cache%h0_Mpc*adotoa)**2*a**5*( Omega_tot_primeprimeprime&
                    & +(Omega_tot*Omega_tot_primeprimeprime +Omega_tot_primeprime*Omega_tot_prime)/(Omega_tot**2 +4._dl*Omega_phi0)**(0.5) +( 8._dl*Omega_tot_prime*Omega_tot_primeprime*Omega_phi0 &
                    &-Omega_tot**2*Omega_tot_prime*Omega_tot_primeprime )/(Omega_tot**2 +4._dl*Omega_phi0)**(1.5) -12._dl*Omega_phi0*Omega_tot*Omega_tot_prime**3/( Omega_tot**2 +4._dl*Omega_phi0 )**(2.5) )

        Omega_tot_primeprimeprimeprime = 360._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-7) +840._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-8)&
                    & +(120._dl*(rhonu_tot+presnu_tot) + (-38._dl/adotoa -9._dl*Hdot/adotoa**3 +Hdotdot/adotoa**4 &
                    & -3._dl*Hdot**2/adotoa**5 )*presnudot_tot +presnudotdot_tot*( 9._dl/adotoa**2 +3._dl*Hdot/adotoa**4 )&
                    & -presnudotdotdot_tot/adotoa**3)/(eft_par_cache%h0_Mpc**2*a**6)

        Hdddd = 14._dl*adotoa*Hddd -71._dl*adotoa**2*Hdotdot +Hdotdot**2/adotoa +3._dl*Hdot*Hddd/adotoa&
                    &+154._dl*adotoa**3*Hdot -14._dl*Hdot*Hdotdot -3._dl*Hdot**2*Hdotdot/adotoa**2 -60._dl*adotoa**5 &
                    &+(eft_par_cache%h0_Mpc**2*adotoa**3*a**6)/(4._dl)*(Omega_tot_primeprimeprimeprime + (Omega_tot_prime*Omega_tot_primeprimeprime +Omega_tot*Omega_tot_primeprimeprimeprime &
                    &+ Omega_tot_primeprimeprime*Omega_tot_prime+(Omega_tot_primeprime)**2)/((Omega_tot**2 +4._dl*Omega_phi0)**(0.5)) -((Omega_tot_primeprimeprime*Omega_tot&
                    &+ Omega_tot_primeprime*Omega_tot_prime)*Omega_tot_prime*Omega_tot)/((Omega_tot**2 +4._dl*Omega_phi0)**(1.5)) +(8._dl*Omega_phi0*(Omega_tot_primeprime)**2 &
                    &+8._dl*Omega_phi0*Omega_tot_prime*Omega_tot_primeprimeprime -2._dl*Omega_tot*Omega_tot_prime**2*Omega_tot_primeprime -Omega_tot**2*( Omega_tot_primeprime)**2 &
                    &-Omega_tot**2*Omega_tot_prime*Omega_tot_primeprimeprime)/((Omega_tot**2 +  4._dl*Omega_phi0)**(1.5))-3._dl*Omega_tot*Omega_tot_prime*(8._dl*Omega_phi0*Omega_tot_prime*Omega_tot_primeprime &
                    &-Omega_tot**2*Omega_tot_prime*Omega_tot_primeprime)/((Omega_tot**2 +  4._dl*Omega_phi0)**(2.5))-12._dl*Omega_phi0*((Omega_tot_prime)**4+3._dl*Omega_tot*(Omega_tot_prime)**2*Omega_tot_primeprime)/&
                    &((Omega_tot**2 +  4._dl*Omega_phi0)**(2.5)) +60._dl*Omega_phi0*(Omega_tot**2*(Omega_tot_prime)**4)/((Omega_tot**2 +4._dl*Omega_phi0)**(3.5)))

          ! compute the psi field and its derivatives
          psi           = self%csi*(eft_par_cache%h0_Mpc)**2*a2/adotoa**2
          psiprime      = 2._dl*self%csi*(eft_par_cache%h0_Mpc)**2*a/adotoa**4*( adotoa**2 -Hdot )
          psiprimeprime = 2._dl*self%csi*(eft_par_cache%h0_Mpc)**2/adotoa**4*( adotoa**2 -3._dl*Hdot &
                    & +4._dl*( Hdot/adotoa )**2 -Hdotdot/adotoa )
          psiprimeprimeprime = -4._dl*Hdot*psiprimeprime/(a*adotoa**2)+ 2._dl*self%csi*eft_par_cache%h0_Mpc**2/(a*adotoa**5)&
                    & *( 2._dl*adotoa*Hdot -3._dl*Hdotdot +9._dl*Hdot*Hdotdot/adotoa**2 &
                    & -8._dl*Hdot**3/adotoa**3-Hddd/adotoa )
          psiprimeprimeprimeprime = -4._dl/(a2*adotoa**4)*psiprimeprime*(Hdotdot*adotoa +3._dl*Hdot**2)&
                    &- psiprimeprimeprime/a*(1._dl +9._dl*Hdot/adotoa**2) +2._dl*self%csi*eft_par_cache%h0_Mpc**2/(a2*adotoa**6)*(2._dl*Hdot**2 &
                    &+2._dl*adotoa*Hdotdot -3._dl*Hddd+ 9._dl*Hdotdot**2/adotoa**2 -42._dl*Hdot**2*Hdotdot/adotoa**3 &
                    &+24._dl*(Hdot/adotoa)**4 +10._dl*Hdot*Hddd/adotoa**2 -Hdddd/adotoa)
          !
          ! scalar field conversion
          phip1 = psi/a
          phip2 = ( psiprime*a -psi )/a2
          phip3 = ( psiprimeprime*a2 -2._dl*psiprime*a +2._dl*psi )/a**3
          phip4 = ( psiprimeprimeprime -3._dl*phip3 )/a
          phip5 = -( psiprimeprimeprime -3._dl*phip3 )/a2+( psiprimeprimeprimeprime -3._dl*phip4 )/a
          !
         ! compute the background EFT functions:
          self%EFTOmega%y(ind)    =-(adotoa**4*phip1**4*(-(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0) &
              &+6*self%QuinticGalileon_c5*(Hdot*phip1 + a*adotoa**2*phip2)))/(2.*a*eft_par_cache%h0_Mpc**6*m0**5)
          !
          self%EFTOmega%yp(ind)    = -((adotoa**2*phip1**3*(-2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(Hdot*phip1 &
                    &+ a*adotoa**2*phip2) +3*self%QuinticGalileon_c5*(4*Hdot**2*phip1**2 + adotoa*phip1**2*Hdotdot &
                    &-adotoa**2*Hdot*phip1*(phip1 - 11*a*phip2) +a**2*adotoa**4*(4*phip2**2 + phip1*phip3))))/(a**2*eft_par_cache%h0_Mpc**6*m0**5))
                    !
          self%EFTOmega%ypp(ind)   = -((phip1**2*(4*Hdot**2*phip1**2*(-(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0) &
                    &+ 6*self%QuinticGalileon_c5*Hdot*phip1) +adotoa*phip1**2*(-2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0 &
                    &+33*self%QuinticGalileon_c5*Hdot*phip1)*Hdotdot +3*adotoa**3*self%QuinticGalileon_c5*phip1**2*Hdotdot*(-3*phip1 + 16*a*phip2) &
                    &+adotoa**2*phip1*(2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdot*(phip1 - 8*a*phip2) &
                    &+3*self%QuinticGalileon_c5*phip1*(Hdot**2*(-12*phip1 + 64*a*phip2) +phip1*Hddd)) +adotoa**4*(&
                    &-2*a**3*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(3*phip2**2 + phip1*phip3) +3*self%QuinticGalileon_c5*Hdot*phip1*(2*phip1**2 &
                    &+ 68*a**2*phip2**2 +a*phip1*(-16*phip2 + 17*a*phip3))) +3*a**3*adotoa**6*self%QuinticGalileon_c5*(12*phip2**3 &
                    &+ 12*phip1*phip2*phip3 +phip1**2*phip4)))/(a**3*eft_par_cache%h0_Mpc**6*m0**5))
                    !
          self%EFTOmega%yppp(ind)  = -((phip1*(5*Hdot*phip1**3*(-2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0 &
                    &+ 21*self%QuinticGalileon_c5*Hdot*phip1)*Hdotdot + adotoa*phip1**2*(-2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(&
                    &-6*Hdot**2*(phip1 - 4*a*phip2) + phip1*Hddd) +3*self%QuinticGalileon_c5*phip1*(11*phip1*Hdotdot**2 &
                    &-24*Hdot**3*(2*phip1 - 7*a*phip2) +13*Hdot*phip1*Hddd)) +3*adotoa**4*self%QuinticGalileon_c5*phip1**2*Hdotdot*(11*phip1**2 &
                    &+ 132*a**2*phip2**2 +3*a*phip1*(-21*phip2 + 11*a*phip3)) +adotoa**3*phip1*(-4*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdot*(phip1**2&
                    & + 18*a**2*phip2**2 +6*a*phip1*(-phip2 + a*phip3)) +3*self%QuinticGalileon_c5*phip1*(3*phip1*(-2*phip1 + 7*a*phip2)*Hddd &
                    &+ 4*Hdot**2*(11*phip1**2 + 132*a**2*phip2**2 +3*a*phip1*(-21*phip2 + 11*a*phip3)))) &
                    &+3*adotoa**2*phip1**2*(2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdotdot*(phip1 - 4*a*phip2) &
                    &+self%QuinticGalileon_c5*phip1*(-33*Hdot*Hdotdot*(2*phip1 - 7*a*phip2) + phip1*Hdddd)) &
                    &-adotoa**5*(2*a**4*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(6*phip2**3 + 9*phip1*phip2*phip3 +phip1**2*phip4) &
                    &+3*self%QuinticGalileon_c5*Hdot*phip1*(6*phip1**3 - 276*a**3*phip2**3 +12*a**2*phip1*phip2*(11*phip2 - 23*a*phip3) +a*phip1**2*(-42*phip2 &
                    &+ 33*a*phip3 - 23*a**2*phip4))) +3*a**4*adotoa**7*self%QuinticGalileon_c5*(24*phip2**4 +72*phip1*phip2**2*phip3 &
                    &+16*phip1**2*phip2*phip4 +phip1**2*(12*phip3**2 + phip1*phip5)) ))/(a**4*adotoa*eft_par_cache%h0_Mpc**6*m0**5))
                    !
          self%EFTc%y(ind)         = (adotoa**2*phip1**2*(a**3*self%QuinticGalileon_c2*eft_par_cache%h0_Mpc**6*m0**3 &
                    &+2*Hdot*phip1*(-(a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2) +3*Hdot*phip1*(&
                    &-(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0) + 6*self%QuinticGalileon_c5*Hdot*phip1)) &
                    &+2*adotoa*phip1**2*(-(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0) &
                    &+18*self%QuinticGalileon_c5*Hdot*phip1)*Hdotdot -12*adotoa**3*self%QuinticGalileon_c5*phip1**2*Hdotdot*(phip1 &
                    &- 4*a*phip2) +adotoa**2*(-2*a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2*(-3*phip1 + a*phip2) &
                    &-3*phip1*(2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdot*(2*phip1 + 3*a*phip2) &
                    &+self%QuinticGalileon_c5*phip1*(15*Hdot**2*(phip1 - 5*a*phip2) -phip1*Hddd))) &
                    &-adotoa**4*(3*self%QuinticGalileon_c5*Hdot*phip1*(6*phip1**2 - 72*a**2*phip2**2 +a*phip1*(25*phip2 &
                    &- 18*a*phip3)) +2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(-11*phip1**2 + 3*a**2*phip2**2 +a*phip1*(5*phip2 &
                    &+ a*phip3))) +3*adotoa**6*self%QuinticGalileon_c5*(7*phip1**3 + 12*a**3*phip2**3 +4*a**2*phip1*phip2*(-phip2 &
                    &+ 3*a*phip3) +a*phip1**2*(-7*phip2 - a*phip3 + a**2*phip4))))/(2.*a*eft_par_cache%h0_Mpc**6*m0**5)

          self%EFTLambda%y(ind)    = (adotoa**2*phip1**2*(a**3*self%QuinticGalileon_c2*eft_par_cache%h0_Mpc**6*m0**3 &
                    &+2*(6*Hdot**2*phip1**2*(-(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0) &
                    &+ 6*self%QuinticGalileon_c5*Hdot*phip1) +2*adotoa*phip1**2*(-(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0) &
                    &+18*self%QuinticGalileon_c5*Hdot*phip1)*Hdotdot -3*adotoa**3*self%QuinticGalileon_c5*phip1**2*Hdotdot*(phip1 &
                    &- 16*a*phip2) -2*a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2*(Hdot*phip1 + a*adotoa**2*phip2) &
                    &+3*adotoa**2*phip1*(-6*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdot*(phip1 + a*phip2) &
                    &+self%QuinticGalileon_c5*phip1*(-3*Hdot**2*(phip1 - 25*a*phip2) +phip1*Hddd)) -2*adotoa**4*(&
                    &a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(phip1**2 + 3*a**2*phip2**2 +a*phip1*(8*phip2 + a*phip3)) &
                    &-3*self%QuinticGalileon_c5*Hdot*phip1*(-3*phip1**2 + 36*a**2*phip2**2 +a*phip1*(4*phip2 + 9*a*phip3))) &
                    &+3*a*adotoa**6*self%QuinticGalileon_c5*(12*a**2*phip2**3 +4*a*phip1*phip2*(2*phip2 + 3*a*phip3) +phip1**2*(&
                    &-4*phip2 + 2*a*phip3 + a**2*phip4)))))/(2.*a*eft_par_cache%h0_Mpc**6*m0**5)

          self%EFTc%yp(ind)      = -(adotoa*phip1*(4*Hdot**2*phip1**2*(a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2 &
                    &+3*Hdot*phip1*(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0 - 6*self%QuinticGalileon_c5*Hdot*phip1)) &
                    &+2*adotoa*phip1**2*(a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2 +9*Hdot*phip1*(&
                    &a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0 - 12*self%QuinticGalileon_c5*Hdot*phip1))*Hdotdot&
                    & - 2*a**3*self%QuinticGalileon_c2*eft_par_cache%h0_Mpc**6*m0**3*(Hdot*phip1 + a*adotoa**2*phip2) &
                    &+2*adotoa**2*phip1*(a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2*Hdot*(-13*phip1 &
                    &+ 7*a*phip2) +phip1*(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(6*Hdot**2*(3*phip1 + 8*a*phip2) &
                    &+ phip1*Hddd)- 6*self%QuinticGalileon_c5*phip1*(3*phip1*Hdotdot**2 +Hdot**3*(-24*phip1 &
                    &+ 90*a*phip2) +4*Hdot*phip1*Hddd))) -3*adotoa**5*self%QuinticGalileon_c5*phip1**2*Hdotdot*(6*phip1**2 &
                    &+ 136*a**2*phip2**2 +a*phip1*(-77*phip2 + 34*a*phip3)) +adotoa**4*(2*a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2*(3*phip1**2 &
                    &+ 2*a**2*phip2**2 +a*phip1*(-9*phip2 + a*phip3)) +3*phip1*(2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdot*(&
                    &-26*phip1**2 + 15*a**2*phip2**2 +5*a*phip1*(3*phip2 + a*phip3)) -self%QuinticGalileon_c5*phip1*(-7*phip1*(phip1 - 3*a*phip2)*Hddd &
                    &+ 3*Hdot**2*(3*phip1**2 + 244*a**2*phip2**2 +a*phip1*(-125*phip2 + 61*a*phip3))))) +adotoa**3*phip1**2*(&
                    &2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdotdot*(4*phip1 + 13*a*phip2) +3*self%QuinticGalileon_c5*phip1*(&
                    &2*Hdot*Hdotdot*(43*phip1 - 145*a*phip2) - phip1*Hdddd)) +adotoa**6*(2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(&
                    &22*phip1**3 + 6*a**3*phip2**3 +3*a**2*phip1*phip2*(5*phip2 + 3*a*phip3) +a*phip1**2*(-49*phip2 + 5*a*phip3 + a**2*phip4)) &
                    &-3*self%QuinticGalileon_c5*Hdot*phip1*(74*phip1**3 + 312*a**3*phip2**3 +12*a**2*phip1*phip2*(-17*phip2 &
                    &+ 26*a*phip3) +a*phip1**2*(-36*phip2 - 51*a*phip3 + 26*a**2*phip4))) +3*adotoa**8*self%QuinticGalileon_c5*(&
                    &21*phip1**4 - 24*a**4*phip2**4 +12*a**3*phip1*phip2**2*(phip2 - 6*a*phip3) -4*a**2*phip1**2*(-6*phip2**2 &
                    &+ 3*a**2*phip3**2 +a*phip2*(-3*phip3 + 4*a*phip4)) +a*phip1**3*(-49*phip2 +a*(6*phip3 + a*(phip4 - a*phip5))))))/(2.*a*eft_par_cache%h0_Mpc**6*m0**5)

          c3= -1.0_dl*self%QuinticGalileon_c3/(eft_par_cache%h0_Mpc)**2

          self%EFTLambda%yp(ind) = self%QuinticGalileon_c2*a2*( adotoa*Hdot*phip1**2 +a*adotoa**3*phip1*phip2 ) &
                    & +4._dl*c3*a2*( adotoa*Hdot*phip1**2 +a*adotoa**3*phip1*phip2 )*( Hdot*phip1/a &
                    & +adotoa**2*phip2 ) +2._dl*c3*a2*adotoa**2*phip1**2*( Hdotdot*phip1/a +adotoa*Hdot*( 3._dl*phip2 &
                    & -phip1/a ) +a*adotoa**3*phip3 ) +1._dl*(+(adotoa*phip1*(4._dl*Hdot**2*phip1**2*( &
                    &+3._dl*Hdot*phip1*(-(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0) + 6._dl*self%QuinticGalileon_c5*Hdot*phip1)) &
                    &-2._dl*adotoa*phip1**2*( +9._dl*Hdot*phip1*(&
                    &a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0 - 12._dl*self%QuinticGalileon_c5*Hdot*phip1))*Hdotdot &
                    &-2._dl*adotoa**2*phip1*( &
                    &+phip1*(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(6._dl*Hdot**2*(5._dl*phip1 + 8._dl*a*phip2) + phip1*Hddd)&
                    &-6._dl*self%QuinticGalileon_c5*phip1*(3._dl*phip1*Hdotdot**2 +Hdot**3*(-12._dl*phip1 +90._dl*a*phip2) &
                    &+4._dl*Hdot*phip1*Hddd))) +3._dl*adotoa**5*self%QuinticGalileon_c5*phip1**2*Hdotdot*(-3._dl*phip1**2 +136._dl*a**2*phip2**2 &
                    &+a*phip1*(-29._dl*phip2 +34._dl*a*phip3)) -adotoa**4*(&
                    & +6._dl*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdot*phip1*(-4._dl*phip1**2 +15._dl*a**2*phip2**2 &
                    &+5._dl*a*phip1*(5._dl*phip2 + a*phip3)) +3._dl*self%QuinticGalileon_c5*phip1**2*(phip1*(4._dl*phip1 -21._dl*a*phip2)*Hddd &
                    &+3._dl*Hdot**2*(9._dl*phip1**2 -244._dl*a**2*phip2**2 +a*phip1*(39._dl*phip2 -61._dl*a*phip3)))) +adotoa**3*phip1**2*(&
                    &-2._dl*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdotdot*(7._dl*phip1 +13._dl*a*phip2) +3._dl*self%QuinticGalileon_c5*phip1*(&
                    &Hdot*Hdotdot*(-47._dl*phip1 +290._dl*a*phip2) + phip1*Hdddd)) -2._dl*adotoa**6*(&
                    &a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(-2._dl*phip1**3 +6._dl*a**3*phip2**3 +3._dl*a**2*phip1*phip2*(8._dl*phip2 +3._dl*a*phip3) &
                    &+a*phip1**2*(-4._dl*phip2 +8._dl*a*phip3 + a**2*phip4)) -3._dl*self%QuinticGalileon_c5*Hdot*phip1*(9._dl*phip1**3 +156._dl*a**3*phip2**3 &
                    &+12._dl*a**2*phip1*phip2*(phip2 +13._dl*a*phip3) +a*phip1**2*(-39._dl*phip2 +3._dl*a*phip3 +13._dl*a**2*phip4))) +3._dl*a*adotoa**8*self%QuinticGalileon_c5*(&
                    &24._dl*a**3*phip2**4 +24._dl*a**2*phip1*phip2**2*(phip2 +3._dl*a*phip3) +4._dl*a*phip1**2*(-6._dl*phip2**2 &
                    &+ 3._dl*a**2*phip3**2 +2._dl*a*phip2*(3._dl*phip3 + 2._dl*a*phip4)) &
                    &+phip1**3*(8._dl*phip2 +a*(-6._dl*phip3 +2._dl*a*phip4 +a**2*phip5)))))/(a*eft_par_cache%h0_Mpc**6*m0**5))

          self%EFTgamma1%y(ind)  = (adotoa**2*phip1**2*(2*Hdot*phip1*(a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2 &
                        &+3*Hdot*phip1*(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0 - 6*self%QuinticGalileon_c5*Hdot*phip1)) &
                        &+2*adotoa*phip1**2*(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0 -18*self%QuinticGalileon_c5*Hdot*phip1)*&
                        &Hdotdot +12*adotoa**3*self%QuinticGalileon_c5*phip1**2*Hdotdot*(phip1 - 4*a*phip2) &
                        &+adotoa**2*(2*a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2*(3*phip1 + a*phip2) +3*phip1*(&
                        &2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdot*(2*phip1 + 3*a*phip2) +self%QuinticGalileon_c5*phip1*(&
                        &15*Hdot**2*(phip1 - 5*a*phip2) -phip1*Hddd))) +adotoa**4*(3*self%QuinticGalileon_c5*Hdot*phip1*(6*phip1**2 &
                        &- 72*a**2*phip2**2 +a*phip1*(25*phip2 - 18*a*phip3)) +2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(16*phip1**2 &
                        &+ 3*a**2*phip2**2 +a*phip1*(5*phip2 + a*phip3))) +3*adotoa**6*self%QuinticGalileon_c5*(13*phip1**3 - 12*a**3*phip2**3 &
                        &+4*a**2*phip1*phip2*(phip2 - 3*a*phip3) +a*phip1**2*(7*phip2 + a*(phip3 - a*phip4)))))/(4.*a**3*eft_par_cache%h0_Mpc**8*m0**5)
                        !
          self%EFTgamma1%yp(ind)  = (phip1*(4*Hdot**2*phip1**2*(a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2 &
                        &+3*Hdot*phip1*(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0 - 6*self%QuinticGalileon_c5*Hdot*phip1)) &
                        &+2*adotoa*phip1**2*(a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2 +9*Hdot*phip1*(&
                        &a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0 - 12*self%QuinticGalileon_c5*Hdot*phip1))*Hdotdot &
                        &+ 2*adotoa**2*phip1*(a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2*Hdot*(11*phip1 &
                        &+ 7*a*phip2) +phip1*(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(6*Hdot**2*(3*phip1 + 8*a*phip2) &
                        &+ phip1*Hddd)- 6*self%QuinticGalileon_c5*phip1*(3*phip1*Hdotdot**2 +Hdot**3*(-24*phip1 &
                        &+ 90*a*phip2) +4*Hdot*phip1*Hddd))) -3*adotoa**5*self%QuinticGalileon_c5*phip1**2*Hdotdot*(&
                        &6*phip1**2 + 136*a**2*phip2**2 +a*phip1*(-77*phip2 + 34*a*phip3)) +adotoa**4*(2*a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2*(&
                        &-3*phip1**2 + 2*a**2*phip2**2 +a*phip1*(9*phip2 + a*phip3)) +3*phip1*(2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdot*(&
                        &28*phip1**2 + 15*a**2*phip2**2 +5*a*phip1*(3*phip2 + a*phip3)) -self%QuinticGalileon_c5*phip1*(-7*phip1*(phip1 &
                        &- 3*a*phip2)*Hddd + 3*Hdot**2*(3*phip1**2 + 244*a**2*phip2**2 +a*phip1*(-125*phip2 + 61*a*phip3))))) &
                        &+adotoa**3*phip1**2*(2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdotdot*(4*phip1 &
                        &+ 13*a*phip2) +3*self%QuinticGalileon_c5*phip1*(2*Hdot*Hdotdot*(43*phip1 - 145*a*phip2) &
                        &- phip1*Hdddd)) +adotoa**6*(3*self%QuinticGalileon_c5*Hdot*phip1*(86*phip1**3 - 312*a**3*phip2**3 &
                        &+12*a**2*phip1*phip2*(17*phip2 - 26*a*phip3) +a*phip1**2*(36*phip2 + 51*a*phip3 - 26*a**2*phip4)) &
                        &+2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(-32*phip1**3 + 6*a**3*phip2**3 +3*a**2*phip1*phip2*(5*phip2 &
                        &+ 3*a*phip3) +a*phip1**2*(59*phip2 + 5*a*phip3 + a**2*phip4))) -3*adotoa**8*self%QuinticGalileon_c5*(39*phip1**4 &
                        &+ 24*a**4*phip2**4 +12*a**3*phip1*phip2**2*(-phip2 + 6*a*phip3) +4*a**2*phip1**2*(-6*phip2**2 + 3*a**2*phip3**2 +a*phip2*(&
                        &-3*phip3 + 4*a*phip4)) +a*phip1**3*(-51*phip2 +a*(-6*phip3 - a*phip4 + a**2*phip5)))))/(4.*a**4*eft_par_cache%h0_Mpc**8*m0**5)

          self%EFTgamma2%y(ind)  = (adotoa**3*phip1**3*(-2._dl*a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2 &
                        &+2._dl*Hdot*phip1*(-(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0) +6._dl*self%QuinticGalileon_c5*Hdot*phip1) &
                        &+3._dl*adotoa*self%QuinticGalileon_c5*phip1**2*Hdotdot -adotoa**2*(3._dl*self%QuinticGalileon_c5*Hdot*phip1*(phip1 &
                        &-11._dl*a*phip2) +2._dl*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(6._dl*phip1 + a*phip2)) +3._dl*adotoa**4*self%QuinticGalileon_c5*(&
                        &-5._dl*phip1**2 +4._dl*a**2*phip2**2 +a**2*phip1*phip3)))/(a**2*eft_par_cache%h0_Mpc**7*m0**5)

          self%EFTgamma2%yp(ind)  = (adotoa*phip1**2*(6._dl*Hdot**2*phip1**2*(-(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0) &
                        &+ 6._dl*self%QuinticGalileon_c5*Hdot*phip1) +2._dl*adotoa*phip1**2*(-(a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0) &
                        &+18._dl*self%QuinticGalileon_c5*Hdot*phip1)*Hdotdot +3._dl*adotoa**3*self%QuinticGalileon_c5*phip1**2*Hdotdot*(-3._dl*phip1 &
                        &+ 16._dl*a*phip2) -6._dl*a**2*self%QuinticGalileon_c3*eft_par_cache%h0_Mpc**4*m0**2*(Hdot*phip1 + a*adotoa**2*phip2) &
                        &+adotoa**2*phip1*(-2._dl*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*Hdot*(29._dl*phip1 + 9._dl*a*phip2) &
                        &+3._dl*self%QuinticGalileon_c5*phip1*(Hdot**2*(-13._dl*phip1 + 75._dl*a*phip2) +phip1*Hddd)) &
                        &-adotoa**4*(2._dl*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(-6._dl*phip1**2 + 3._dl*a**2*phip2**2 +a*phip1*(24._dl*phip2 &
                        &+ a*phip3)) +3._dl*self%QuinticGalileon_c5*Hdot*phip1*(33._dl*phip1**2 - 72._dl*a**2*phip2**2 -2._dl*a*phip1*(-8._dl*phip2 + 9._dl*a*phip3))) &
                        &+3._dl*adotoa**6*self%QuinticGalileon_c5*(10._dl*phip1**3 +12._dl*a**3*phip2**3 +12._dl*a**3*phip1*phip2*phip3 +phip1**2*(-25._dl*a*phip2 &
                        &+ a**3*phip4))))/(a**3*eft_par_cache%h0_Mpc**7*m0**5)

          self%EFTgamma3%y(ind)  = (adotoa**4*phip1**4*(-2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0 &
                        &+3*self%QuinticGalileon_c5*(Hdot*phip1 + adotoa**2*(-phip1 + a*phip2))))/(a*eft_par_cache%h0_Mpc**6*m0**5)

          self%EFTgamma3%yp(ind)  = (adotoa**2*phip1**3*(-8*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(Hdot*phip1 &
                        &+ a*adotoa**2*phip2) +3*self%QuinticGalileon_c5*(4*Hdot**2*phip1**2 + adotoa*phip1**2*Hdotdot &
                        &+adotoa**2*Hdot*phip1*(-7*phip1 + 11*a*phip2) +adotoa**4*(phip1**2 + 4*a**2*phip2**2 &
                        &+a*phip1*(-5*phip2 + a*phip3)))))/(a**2*eft_par_cache%h0_Mpc**6*m0**5)

          self%EFTgamma4%y(ind)  = -self%EFTgamma3%y(ind)
          self%EFTgamma4%yp(ind)  = -self%EFTgamma3%yp(ind)

          self%EFTgamma4%ypp(ind) =-(phip1**2*(8*Hdot**2*phip1**2*(-2*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0 &
                        &+ 3*self%QuinticGalileon_c5*Hdot*phip1) +adotoa*phip1**2*(-8*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0 &
                        &+33*self%QuinticGalileon_c5*Hdot*phip1)*Hdotdot +3*adotoa**3*self%QuinticGalileon_c5*phip1**2*&
                        &Hdotdot*(-9*phip1 + 16*a*phip2) +adotoa**2*phip1*(8*a*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*&
                        &m0*Hdot*(phip1 - 8*a*phip2) +3*self%QuinticGalileon_c5*phip1*(Hdot**2*(-36*phip1 + 64*a*phip2) &
                        &+phip1*Hddd)) +adotoa**4*(-8*a**3*self%QuinticGalileon_c4*eft_par_cache%h0_Mpc**2*m0*(3*phip2**2 &
                        &+ phip1*phip3) +3*self%QuinticGalileon_c5*Hdot*phip1*(20*phip1**2 + 68*a**2*phip2**2 +a*phip1*(-76*phip2 &
                        &+ 17*a*phip3))) +3*adotoa**6*self%QuinticGalileon_c5*(-2*phip1**3 + 12*a**3*phip2**3 +4*a**2*phip1*phip2*(&
                        &-5*phip2 + 3*a*phip3) +a*phip1**2*(10*phip2 - 5*a*phip3 + a**2*phip4))))/(a**3*eft_par_cache%h0_Mpc**6*m0**5)

          end subroutine

    end subroutine EFTCAMBQuinticGalileonSolveBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBQuinticGalileonComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_Quintic_Galileon)  :: self   !< the base class

        self%parameter_number = 2

    end subroutine EFTCAMBQuinticGalileonComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBQuinticGalileonFeedback( self, print_params )

        implicit none

        class(EFTCAMB_Quintic_Galileon)    :: self         !< the base class
        logical, optional                  :: print_params !< optional flag that decised whether to print numerical values
                                                       !! of the parameters.

        logical                            :: print_params_temp

        ! print general model informations:
        if (self%QuinticGalileon_c2 == -1._dl) then

          write(*,*)
          write(*,'(a,a)')    '   Model               =  ', self%name
          write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
          write(*,'(a,F12.6)')   '                 xi    ='  , self%csi
          write(*,'(a,F12.6)')   '                 c2    ='  , self%QuinticGalileon_c2
          write(*,'(a,F12.6)')   '                 c3    ='  , self%QuinticGalileon_c3
          write(*,'(a,F12.6)')   '                 c4    ='  , self%QuinticGalileon_c4
          write(*,'(a,F12.6)')   '                 c5    ='  , self%QuinticGalileon_c5

        end if
        ! print the values of the parameters:
        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

    end subroutine EFTCAMBQuinticGalileonFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBQuinticGalileonParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_Quintic_Galileon)   :: self   !< the base class
        integer     , intent(in)          :: i      !< the index of the parameter
        character(*), intent(out)         :: name   !< the output name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            name = 'Quintic_Galileon_xi'
            return
        end if
        if ( i==2 ) then
            name = 'Quintic_Galileon_c3'
            return
        end if
        if ( i==0 ) then
            name = 'no_name'
            return
        end if

    end subroutine EFTCAMBQuinticGalileonParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBQuinticGalileonParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_Quintic_Galileon)   :: self       !< the base class
        integer     , intent(in)          :: i          !< The index of the parameter
        character(*), intent(out)         :: latexname  !< the output latex name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names_latex.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            latexname = '\xi'
            return
        end if
        if ( i==2 ) then
            latexname = 'c_3'
            return
        end if
        if ( i==0 ) then
            latexname = 'noname'
            return
        end if

    end subroutine EFTCAMBQuinticGalileonParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBQuinticGalileonParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_Quintic_Galileon)   :: self   !< the base class
        integer , intent(in)              :: i      !< The index of the parameter
        real(dl), intent(out)             :: value  !< the output value of the i-th parameter

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
        if ( i==1 ) then
            value = self%csi
            return
        end if
        if ( i==2 ) then
            value = self%QuinticGalileon_c3
            return
        end if

    end subroutine EFTCAMBQuinticGalileonParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBQuinticGalileonBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Quintic_Galileon)              :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.
        real(dl) :: x, mu
        integer  :: ind

        x   = log(a)
        if(x>=self%x_final) return
        if(x<=self%x_initial) return

        if(a==0._dl) then
            return
        else if (eft_cache%adotoa==0._dl) then
          ! call abort()
            call self%compute_adotoa( a, eft_par_cache, eft_cache )
            call self%compute_H_derivs( a, eft_par_cache, eft_cache )
            if(eft_cache%adotoa==0._dl) return
            if(eft_cache%Hdot==0._dl) return
            if(eft_cache%Hdotdot==0._dl) return
        end if

        call self%EFTOmega%precompute(x, ind, mu )

        ! compute the background EFT functions:
        eft_cache%EFTOmegaV    = self%EFTOmega%value( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaP    = self%EFTOmega%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaPP   = self%EFTOmega%second_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaPPP  = self%EFTOmega%third_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTc         = self%EFTc%value( x, index=ind, coeff=mu )
        eft_cache%EFTLambda    = self%EFTLambda%value( x, index=ind, coeff=mu )
        eft_cache%EFTcdot      = self%EFTc%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTLambdadot = self%EFTLambda%first_derivative( x, index=ind, coeff=mu )

    end subroutine EFTCAMBQuinticGalileonBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the second order EFT functions at a given time.
    subroutine EFTCAMBQuinticGalileonSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Quintic_Galileon)              :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu
        integer  :: ind

        x   = log(a)
        if(x>=self%x_final) return
        if(x<=self%x_initial) return

        call self%EFTgamma1%precompute(x, ind, mu )
        !
        ! ! compute the second order EFT functions:
        eft_cache%EFTGamma1V  = self%EFTgamma1%value( x, index=ind, coeff=mu )
        eft_cache%EFTGamma1P  = self%EFTgamma1%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTGamma2V  = self%EFTgamma2%value( x, index=ind, coeff=mu )
        eft_cache%EFTGamma2P  = self%EFTgamma2%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTGamma3V  = self%EFTgamma3%value( x, index=ind, coeff=mu )
        eft_cache%EFTGamma3P  = self%EFTgamma3%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTGamma4V  = self%EFTgamma4%value( x, index=ind, coeff=mu )
        eft_cache%EFTGamma4P  = self%EFTgamma4%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTGamma4PP = self%EFTgamma4%second_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTGamma5V  = 0.5_dl*eft_cache%EFTGamma3V
        eft_cache%EFTGamma5P  = 0.5_dl*eft_cache%EFTGamma3P
        eft_cache%EFTGamma6V  = 0._dl
        eft_cache%EFTGamma6P  = 0._dl

    end subroutine EFTCAMBQuinticGalileonSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H.
    subroutine EFTCAMBQuinticGalileonComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Quintic_Galileon)              :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp, a2, Omega_tot

        a2 = a*a
        Omega_tot = ( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-3) + ( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-4) +eft_cache%grhonu_tot/(3._dl*eft_par_cache%h0_Mpc**2*a2)
        temp = 0.5_dl*a2*(eft_par_cache%h0_Mpc)**2*( Omega_tot + sqrt( Omega_tot**2 +4._dl*eft_par_cache%omegav ) )
        eft_cache%adotoa = sqrt( temp )

    end subroutine EFTCAMBQuinticGalileonComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBQuinticGalileonComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Quintic_Galileon)              :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp, a2, Omega_tot, Omega_tot_prime, Omega_tot_primeprime,Omega_tot_primeprimeprime, Omega_phi0,Omega_tot_primeprimeprimeprime
        !
        a2 = a*a
        !
        if(a*eft_cache%adotoa==0._dl) return

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

    end subroutine EFTCAMBQuinticGalileonComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBQuinticGalileonAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Quintic_Galileon)              :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBQuinticGalileonAdditionalModelStability       !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBQuinticGalileonAdditionalModelStability = .True.

    end function EFTCAMBQuinticGalileonAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_full_Quintic_Galileon

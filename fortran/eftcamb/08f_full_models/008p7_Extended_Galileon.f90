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

!> @file 10p5_Extended_Galileon.f90
!! This file contains the definition of the Extended Galileon model.
!! Please refer to the numerical notes for details.

!----------------------------------------------------------------------------------------
!> This module contains the definition of the Extended Galileon model.
!! Please refer to the numerical notes for details.

!> @author Luis Atayde, Noemi Frusciante, Simone Peirone

!> updated to pyEFTCAMB by Fabrizio Renzi

module EFTCAMB_full_Extended_Galileon


    use precision
    use IniObjects
    use MpiUtils
    use FileUtils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_model_full

    implicit none

    private

    public EFTCAMB_Extended_Galileon

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of the Extended Galileon model.
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_Extended_Galileon

        ! the model parameters:
        real(dl)  :: ExtendedGalileon_B      !< Extended Galileon model parameter \f$B\f$
        real(dl)  :: ExtendedGalileon_q      !< Extended Galileon model parameter \f$q\f$
        real(dl)  :: csi                     !< Extended Galileon background parameter \f$\xi\f$ deriving from the tracker solution
	    real(dl)  :: S

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBExtendedGalileonReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBExtendedGalileonAllocateModelSelection      !< subroutine that allocates the model selection. For Horava this is a dummy procedure.
        procedure :: init_model_parameters           => EFTCAMBExtendedGalileonInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBExtendedGalileonInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBExtendedGalileonInitBackground              !< subroutine that initializes the background of Extended Galileon.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBExtendedGalileonComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBExtendedGalileonFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBExtendedGalileonParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBExtendedGalileonParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBExtendedGalileonParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBExtendedGalileonBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBExtendedGalileonSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_adotoa                    => EFTCAMBExtendedGalileonComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBExtendedGalileonComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.

        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBExtendedGalileonAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_Extended_Galileon

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBExtendedGalileonReadModelSelectionFromFile( self, Ini, eft_error)

        implicit none

        class(EFTCAMB_Extended_Galileon)       :: self          !< the base class
        type(TIniFile)                         :: Ini           !< Input ini file
        integer                                :: eft_error     !< error code: 0 all fine, 1 initialization failed

    end subroutine EFTCAMBExtendedGalileonReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection. 
    subroutine EFTCAMBExtendedGalileonAllocateModelSelection( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_Extended_Galileon) :: self        !< the base class
        type(TIniFile)                   :: Ini         !< Input ini file
        integer                          :: eft_error   !< error code: 0 all fine, 1 initialization failed

        ! Nothing needs to be done but procedure present because its definition is deferred.

    end subroutine EFTCAMBExtendedGalileonAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBExtendedGalileonInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_Extended_Galileon)                          :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in)    :: array  !< input array with the values of the parameters of the model.

        self%S                  = array(1)
        self%ExtendedGalileon_q = array(2)

    end subroutine EFTCAMBExtendedGalileonInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file. 
    subroutine EFTCAMBExtendedGalileonInitModelParametersFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_Extended_Galileon)  :: self      !< the base class
        type(TIniFile)                    :: Ini       !< Input ini file
        integer                           :: eft_error !< error code: 0 all fine, 1 initialization failed

        self%S                  = Ini%Read_Double('ExtendedGalileon_S', 0._dl )
        self%ExtendedGalileon_q = Ini%Read_Double('ExtendedGalileon_q', 0._dl )

    end subroutine EFTCAMBExtendedGalileonInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBExtendedGalileonComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_Extended_Galileon)  :: self               !< the base class

        self%parameter_number = 2								!luis-change::0->2

    end subroutine EFTCAMBExtendedGalileonComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBExtendedGalileonFeedback( self, print_params )

        implicit none

        class(EFTCAMB_Extended_Galileon)  :: self         !< the base class
        logical, optional                 :: print_params !< optional flag that decised whether to print numerical values
                                                          !! of the parameters.

        ! print general model informations:
       if (self%ExtendedGalileon_B /= 0._dl) then

        write(*,*)
        write(*,'(a,a)')       '   Model              ='  , self%name
        write(*,'(a,I3)')      '   Number of params   ='  , self%parameter_number
        write(*,'(a,F12.6)')   '                xi    ='  , self%csi
        write(*,'(a,F12.6)')   '                B     ='  , self%ExtendedGalileon_B
        write(*,'(a,F12.6)')   '                q     ='  , self%ExtendedGalileon_q
        write(*,'(a,F12.6)')   '                S     ='  , self%S

       end if

    end subroutine EFTCAMBExtendedGalileonFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBExtendedGalileonParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_Extended_Galileon) :: self   !< the base class
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
            name = 'ExtendedGalileon_S'
            return
        end if
        if ( i==2 ) then
            name = 'ExtendedGalileon_q'
            return
        end if
        if ( i==0 ) then
            name = 'noname'
            return
        end if

    end subroutine EFTCAMBExtendedGalileonParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBExtendedGalileonParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_Extended_Galileon) :: self       !< the base class
        integer     , intent(in)         :: i          !< The index of the parameter
        character(*), intent(out)        :: latexname  !< the output latex name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names_latex.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            latexname = 'S'
            return
        end if
        if ( i==2 ) then
            latexname = 'q'
            return
        end if
        if ( i==0 ) then
            latexname = 'noname'
            return
        end if

    end subroutine EFTCAMBExtendedGalileonParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBExtendedGalileonParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_Extended_Galileon) :: self   !< the base class
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
        if ( i==1 ) then
            value = self%S
            return
        end if
        if ( i==2 ) then
            value = self%ExtendedGalileon_q
            return
        end if

    end subroutine EFTCAMBExtendedGalileonParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of Extended Galileon.
    subroutine EFTCAMBExtendedGalileonInitBackground( self, params_cache, feedback_level, success, outroot  )

        implicit none

        class(EFTCAMB_Extended_Galileon)              :: self           !< the base class
        type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                       , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                       , intent(out)   :: success        !< wether the background initialization succeded or not
        character(LEN=*), optional    , intent(in)    :: outroot        !< the output root for the debug files

        real(dl) :: Omega_phi0

        Omega_phi0 = params_cache%omegav

        self%ExtendedGalileon_B = self%S*self%ExtendedGalileon_q			!< B=S*q!

        call self%feedback()

        success=.true.

    end subroutine EFTCAMBExtendedGalileonInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBExtendedGalileonBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache)

        implicit none

        class(EFTCAMB_Extended_Galileon)                :: self          !< the base class
        real(dl), intent(in)                            :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout)   :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout)   :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.


        real(dl)    :: a2, Omega_phi0

        Omega_phi0 = eft_par_cache%omegav
        a2 = a*a

        if(a==0._dl) then
            return
        else if (eft_cache%adotoa==0._dl) then
            call self%compute_adotoa( a, eft_par_cache, eft_cache )
            call self%compute_H_derivs( a, eft_par_cache, eft_cache )
        end if

        if (eft_cache%adotoa <= 0._dl .or. a <= 0._dl) return
        ! compute the background EFT functions:
        eft_cache%EFTOmegaV    = 0._dl
        eft_cache%EFTOmegaP    = 0._dl
        eft_cache%EFTOmegaPP   = 0._dl
        eft_cache%EFTOmegaPPP  = 0._dl
        eft_cache%EFTc         = -(a**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*self%ExtendedGalileon_B&
	                                 &*eft_cache%adotoa**(-2 - self%ExtendedGalileon_B/self%ExtendedGalileon_q)&
	                                 &*eft_par_cache%h0_Mpc**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)&
	                                 &*Omega_phi0*(eft_cache%adotoa**2 - eft_cache%Hdot))/(2.*self%ExtendedGalileon_q)


        eft_cache%EFTLambda    = -((a**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)&
	                             &*eft_cache%adotoa**(-2 - self%ExtendedGalileon_B/self%ExtendedGalileon_q)&
	                             &*eft_par_cache%h0_Mpc**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*Omega_phi0*((3*self%ExtendedGalileon_q + self%ExtendedGalileon_B)&
	                             &*eft_cache%adotoa**2 - self%ExtendedGalileon_B*eft_cache%Hdot))/self%ExtendedGalileon_q)

        eft_cache%EFTcdot      =-(a**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*self%ExtendedGalileon_B&
	                           &*eft_cache%adotoa**(-3 - self%ExtendedGalileon_B/self%ExtendedGalileon_q)*eft_par_cache%h0_Mpc**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)&
	                           &*Omega_phi0*(self%ExtendedGalileon_B*eft_cache%adotoa**4 - 2*self%ExtendedGalileon_B&
	                           &*eft_cache%adotoa**2*eft_cache%Hdot + (2*self%ExtendedGalileon_q + self%ExtendedGalileon_B)*eft_cache%Hdot**2 - self%ExtendedGalileon_q&
	                           &*eft_cache%adotoa*eft_cache%Hdotdot))/(2.*self%ExtendedGalileon_q**2)

        eft_cache%EFTLambdadot = -((a**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)&
	                            &*self%ExtendedGalileon_B*eft_cache%adotoa**(-3 - self%ExtendedGalileon_B/self%ExtendedGalileon_q)&
	                            &*eft_par_cache%h0_Mpc**(2 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*Omega_phi0&
	                            &*((3*self%ExtendedGalileon_q + self%ExtendedGalileon_B)*eft_cache%adotoa**4 - (3*self%ExtendedGalileon_q + 2*self%ExtendedGalileon_B)*eft_cache%adotoa**2&
	                            &*eft_cache%Hdot + (2*self%ExtendedGalileon_q + self%ExtendedGalileon_B)*eft_cache%Hdot**2 - self%ExtendedGalileon_q*eft_cache%adotoa&
	                            &*eft_cache%Hdotdot))/self%ExtendedGalileon_q**2)

    end subroutine EFTCAMBExtendedGalileonBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the second order EFT functions at a given time.
    subroutine EFTCAMBExtendedGalileonSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Extended_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: a2,Omega_phi0

        Omega_phi0 = eft_par_cache%omegav
        a2 = a*a

        if(a*eft_cache%adotoa==0._dl) return
        if (eft_cache%adotoa <= 0._dl .or. a <= 0._dl) return

        ! compute the second order EFT functions:
        eft_cache%EFTGamma1V  = (a**(self%ExtendedGalileon_B/self%ExtendedGalileon_q)*self%ExtendedGalileon_B*eft_cache%adotoa**(-2 - self%ExtendedGalileon_B/self%ExtendedGalileon_q)&
	                            *eft_par_cache%h0_Mpc**(self%ExtendedGalileon_B/self%ExtendedGalileon_q)*Omega_phi0*((1 + 12*self%ExtendedGalileon_q**2)*eft_cache%adotoa**2 - eft_cache%Hdot))/(4.*self%ExtendedGalileon_q)

        eft_cache%EFTGamma1P  =  (a**(-1 + self%ExtendedGalileon_B/self%ExtendedGalileon_q)*self%ExtendedGalileon_B*eft_cache%adotoa**(-4 - self%ExtendedGalileon_B/self%ExtendedGalileon_q)                    &
                                 *eft_par_cache%h0_Mpc**(self%ExtendedGalileon_B/self%ExtendedGalileon_q)*Omega_phi0*((self%ExtendedGalileon_B + 12*self%ExtendedGalileon_q**2*self%ExtendedGalileon_B)         &
                             	 *eft_cache%adotoa**4 - 2*(1 + 6*self%ExtendedGalileon_q**2)*self%ExtendedGalileon_B*eft_cache%adotoa**2*eft_cache%Hdot + (2*self%ExtendedGalileon_q + self%ExtendedGalileon_B) &
                             	 *eft_cache%Hdot**2 - self%ExtendedGalileon_q*eft_cache%adotoa*eft_cache%Hdotdot))/(4.*self%ExtendedGalileon_q**2)

        eft_cache%EFTGamma2V  = -2*self%ExtendedGalileon_B*(a*eft_par_cache%h0_Mpc/eft_cache%adotoa)**((self%ExtendedGalileon_B+self%ExtendedGalileon_q)/self%ExtendedGalileon_q)*Omega_phi0

        eft_cache%EFTGamma2P  = -2/self%ExtendedGalileon_q*self%ExtendedGalileon_B*(self%ExtendedGalileon_B+self%ExtendedGalileon_q)*a**(self%ExtendedGalileon_B/self%ExtendedGalileon_q)&
		                        *eft_par_cache%h0_Mpc**((self%ExtendedGalileon_B+self%ExtendedGalileon_q)/self%ExtendedGalileon_q)&
                                *eft_cache%adotoa**(-3-(self%ExtendedGalileon_B)/self%ExtendedGalileon_q)*Omega_phi0*(eft_cache%adotoa**2-eft_cache%Hdot)

        eft_cache%EFTGamma3V  = 0._dl
        eft_cache%EFTGamma3P  = 0._dl
        eft_cache%EFTGamma4V  = 0._dl
        eft_cache%EFTGamma4P  = 0._dl
        eft_cache%EFTGamma4PP = 0._dl
        eft_cache%EFTGamma5V  = 0._dl
        eft_cache%EFTGamma5P  = 0._dl
        eft_cache%EFTGamma6V  = 0._dl
        eft_cache%EFTGamma6P  = 0._dl

    end subroutine EFTCAMBExtendedGalileonSecondOrderEFTFunctions


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H.
    subroutine EFTCAMBExtendedGalileonComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Extended_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp, a2, Omega_tot
        integer     :: nu_i , counter
	      real(dl)::limit1, limit2, flimit1, flimit2, dmean, soluction, fsoluction, bolean !soluction=H/H0
        real(dl) :: Omega_phi0
	      real(dl) :: ATemp1, ATemp2, BTemp1, BTemp2, HorizAsyntB

        Omega_phi0 = eft_par_cache%omegav
        a2 = a*a

        Omega_tot = ( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-3.) + ( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-4.) +eft_cache%grhonu_tot/(3._dl*eft_par_cache%h0_Mpc**2.*a2)

        !SP:Debug
        if (a<0.2_dl)then
          ! temp = 0.5_dl*a2*(eft_par_cache%h0_Mpc)**2.*( Omega_tot + sqrt( Omega_tot**2. +4._dl*eft_par_cache%omegav ) )
          ! eft_cache%adotoa =sqrt( temp )
          ! if (eft_cache%grhom_t == 0._dl) call self%compute_background_EFT_functions(a, eft_par_cache, eft_cache )
          temp = 1.0_dl/(1.0_dl + eft_cache%EFTOmegaV+ a*eft_cache%EFTOmegaP)*(eft_cache%grhom_t + 2.0_dl*eft_cache%EFTc -eft_cache%EFTLambda )/3.0_dl
          eft_cache%adotoa = sqrt(temp)
          return
        end if

        limit1=0
        if (limit1.lt.0) limit1=0
          limit2=10**(9)
          flimit1=Omega_phi0+Omega_tot*(limit1/a)**(self%S)-(limit1/a)**(2+self%S)
          flimit2=Omega_phi0+Omega_tot*(limit2/a)**(self%S)-(limit2/a)**(2+self%S)
          dmean=(limit2-limit1)/2
          soluction=limit2-dmean
          fsoluction=1
          counter=0
          do while(sqrt(fsoluction**2).gt.10**(-1).and.counter.lt.50**1)
	             fsoluction=Omega_phi0+Omega_tot*(soluction/a)**(self%S)-(soluction/a)**(2+self%S)
	             bolean=fsoluction*flimit1
	             if (bolean.gt.0.) then
		               limit1=soluction
		               flimit1=fsoluction
	             endif
	             if (bolean.le.0.) then
		               limit2=soluction
		               flimit2=fsoluction
	             endif
	             dmean=(limit2-limit1)/2
	             soluction=limit1+dmean
	             counter=counter+1
          enddo

		  temp= soluction*eft_par_cache%h0_Mpc
      eft_cache%adotoa = temp

    end subroutine EFTCAMBExtendedGalileonComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBExtendedGalileonComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Extended_Galileon)             :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp, a2, Omega_tot, Omega_tot_prime, Omega_tot_primeprime, Omega_phi0
        integer     :: nu_i

        a2 = a*a

        if (a == 0._dl) return
        if (eft_cache%adotoa == 0._dl) then
          call self%compute_adotoa( a, eft_par_cache, eft_cache )
          ! print*, "here = ", eft_cache%adotoa, a
          if (eft_cache%adotoa == 0._dl) return
        end if

        Omega_tot = ( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-3) + ( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-4) +eft_cache%grhonu_tot/(3._dl*eft_par_cache%h0_Mpc**2*a2)

        Omega_tot_prime = -3._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-4) -4._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-5) &
                          -(eft_cache%grhonu_tot+eft_cache%gpinu_tot)/(eft_par_cache%h0_Mpc**2*a2*a)

        Omega_tot_primeprime = 12._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-5) +20._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-6)&
                             +(4._dl*(eft_cache%grhonu_tot+eft_cache%gpinu_tot)-eft_cache%gpinudot_tot/eft_cache%adotoa )/(eft_par_cache%h0_Mpc**2*a2**2)

        Omega_phi0 = eft_par_cache%omegav

        eft_cache%Hdot = -(eft_cache%adotoa**2*(a**2*eft_par_cache%h0_Mpc**2*(a*Omega_tot_prime-self%S*Omega_tot)+(self%S+2)*eft_cache%adotoa**2))/(a**2*eft_par_cache%h0_Mpc**2*self%S*Omega_tot-(self%S+2)*eft_cache%adotoa**2)
        eft_cache%Hdotdot = eft_cache%adotoa**3/((self%S+2)*eft_cache%adotoa**2*a**2*eft_par_cache%h0_Mpc**2*self%S*Omega_tot)**3*(a**2*eft_par_cache%h0_Mpc**2*(self%S+2)**2&
                          * eft_cache%adotoa**4*(a*(5*Omega_tot_prime+a*Omega_tot_primeprime)-6*self%S*Omega_tot)+a**6*eft_par_cache%h0_Mpc**6*self%S*Omega_tot*(-a**2*(self%S+2)*&
                            Omega_tot_prime**2-2*self%S**2*Omega_tot**2+a*self%S*Omega_tot*(5*Omega_tot_prime+a*Omega_tot_primeprime))+a**4*eft_par_cache%h0_Mpc**4*self%S*(self%S+2)*eft_cache%adotoa**2*(a**2*Omega_tot_prime**2&
                          + 6*self%S*Omega_tot**2-2*a*Omega_tot*(5*Omega_tot_prime+a*Omega_tot_primeprime))+2*(self%S+2)**3*eft_cache%adotoa**6)

    end subroutine EFTCAMBExtendedGalileonComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBExtendedGalileonAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Extended_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBExtendedGalileonAdditionalModelStability       !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBExtendedGalileonAdditionalModelStability = .True.

	!mass condition
	if ( (self%S+0.02)*self%ExtendedGalileon_q<.5 ) EFTCAMBExtendedGalileonAdditionalModelStability = .false.

	!Strong Coupling
	if ( (self%ExtendedGalileon_q*self%S -1.)/(2.*self%ExtendedGalileon_q*(1. +self%S/2.))>0.) EFTCAMBExtendedGalileonAdditionalModelStability = .false.

	return

    end function EFTCAMBExtendedGalileonAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_full_Extended_Galileon

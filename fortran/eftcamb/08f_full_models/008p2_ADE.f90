!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2023 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 008p5_ADE.f90
!! This file contains the implementation of the Acoustic Dark Energy model discussed in

!----------------------------------------------------------------------------------------
!> This module contains

!> @author Giampaolo Benevento, Marco Raveri

module EFTCAMB_ADE_Mod

    use precision
    use IniObjects
    use MpiUtils
    use FileUtils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_model_full
    use equispaced_linear_interpolation_1D
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_constant_parametrization_1D
    use EFTCAMB_quadpack
    use MassiveNu

    implicit none

    private

    public EFTCAMB_ADE

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of ADE models.
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_ADE


        ! the model parameters:
        real(dl)  :: cs2                 !< The integer from which the ADE equation of state is computed
        real(dl)  :: Log_ac              !< Reference scale factor for ADE
        real(dl)  :: f_ac                !< Density parameter of ADE at a=ac defined as rho_ADE(a_c)/rho_tot(a_c)
        real(dl)  :: p                   !< ADE exponential parameter
        real(dl)  :: wf                  !< value of w after a_c

        ! the interpolated EFT functions that come out of the background sover:
        type(equispaced_linear_interpolate_function_1D) :: EFTLambda      !< The interpolated function Lambda (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTc           !< The interpolated function c (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTgamma1      !< The interpolated function gamma1 (and derivatives).

        ! some designer parameters:
        integer  :: background_num_points = 1500                          !< Number of points sampled by the designer code.
        real(dl) :: x_initial             = log(10._dl**(-10._dl))        !< log(a start)
        real(dl) :: x_final               = 0.0_dl                        !< log(a final)

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBADEReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBADEAllocateModelSelection      !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBADEInitModelParameters         !< subroutine that initializes the model
        procedure :: init_model_parameters_from_file => EFTCAMBADEInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBADEComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBADEFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBADEParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBADEParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBADEParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBADEBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBADESecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBADEInitBackground             !< subroutine that initializes the background of ADE.
        procedure :: solve_background_equations      => EFTCAMBADESolveBackgroundEquations   !< subroutine that solves the background equations.
        ! stability procedures:
        procedure :: additional_model_stability      => EFTCAMBADEAdditionalModelStability   !< function that computes model specific stability requirements.



    end type EFTCAMB_ADE

    ! ---------------------------------------------------------------------------------------------

    ! define debug files
    type(TTextFile) :: file_debug_1, file_debug_2

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBADEReadModelSelectionFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_ADE) :: self      !< the base class
        type(TIniFile)     :: Ini       !< Input ini file
        integer            :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! read model selection flags:

        ! read precision parameters
        self%background_num_points = Ini%Read_Int( 'model_background_num_points', 1500 )
        self%x_initial = Log( Ini%Read_Double( 'model_background_a_ini', 1d-10 ) )
        self%x_final = Log( Ini%Read_Double( 'model_background_a_final', 1._dl ) )

    end subroutine EFTCAMBADEReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBADEAllocateModelSelection( self, ini, eft_error )

        implicit none

        class(EFTCAMB_ADE) :: self      !< the base class
        type(TIniFile)     :: Ini       !< Input ini file
        integer            :: eft_error !< error code: 0 all fine, 1 initialization failed

        ! nothing to be done here.

    end subroutine EFTCAMBADEAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBADEInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_ADE)                                     :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

        self%cs2     = array(1)
        self%Log_ac  = array(2)
        self%f_ac    = array(3)
        self%p       = array(4)
        self%wf      = array(5)

    end subroutine EFTCAMBADEInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBADEInitModelParametersFromFile( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_ADE) :: self      !< the base class
        type(TIniFile)     :: Ini       !< Input ini file
        integer            :: eft_error !< error code: 0 all fine, 1 initialization failed

        self%cs2    = Ini%Read_Double('cs2'   , 1._dl  )
        self%Log_ac = Ini%Read_Double('Log_ac', -4._dl )
        self%f_ac   = Ini%Read_Double('f_ac'  , 0.07_dl)
        self%p      = Ini%Read_Double('p'     , 1._dl  )
        self%wf     = Ini%Read_Double('wf'    , 1._dl  )

    end subroutine EFTCAMBADEInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBADEComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_ADE)  :: self   !< the base class

        self%parameter_number = 5

    end subroutine EFTCAMBADEComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBADEFeedback( self, print_params )

        implicit none

        class(EFTCAMB_ADE)  :: self         !< the base class
        logical, optional   :: print_params !< optional flag that decised whether to print numerical values
                                            !! of the parameters.

        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name

        ! print model informations:

        write(*,*)
        write(*,'(a24,F12.6)') 'cs2         =', self%cs2
        write(*,'(a24,F12.6)') 'a_c         =', 10**self%Log_ac
        write(*,'(a24,F12.6)') 'f_ADE(a_c)  =', self%f_ac
        write(*,'(a24,F12.6)') 'p           =', self%p
        write(*,'(a24,F12.6)') 'wf          =', self%wf

    end subroutine EFTCAMBADEFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBADEParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_ADE)    :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! other parameters:
        else if ( i==1 ) then
            name = TRIM('cs2')
        else if ( i==2 ) then
            name = TRIM('Log_ac')
        else if ( i==3 ) then
            name = TRIM('f_ac')
        else if ( i==4 ) then
            name = TRIM('p')
        else if ( i==5 ) then
            name = TRIM('wf')
            return
        end if

    end subroutine EFTCAMBADEParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBADEParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_ADE)    :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! other parameters:
        else if ( i==1 ) then
            latexname = TRIM('c_s^2')
            return
        else if ( i==2 ) then
            latexname = TRIM('\\log(a_c)')
            return
        else if ( i==3 ) then
            latexname = TRIM('f(a_c)')
            return
        else if ( i==4 ) then
            latexname = TRIM('p')
            return
        else if ( i==5 ) then
            latexname = TRIM('w_f')
            return
        end if

    end subroutine EFTCAMBADEParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBADEParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_ADE)    :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! other parameters:
        else if ( i==1 ) then
            value = self%cs2
            return
        else if ( i==2 ) then
            value = self%Log_ac
            return
        else if ( i==3 ) then
            value = self%f_ac
            return
        else if ( i==4 ) then
            value = self%p
            return
        else if ( i==5 ) then
            value = self%wf
            return
        end if

    end subroutine EFTCAMBADEParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBADEBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_ADE)                            :: self          !< the base class
        real(dl), intent(in)                          :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout)  :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu
        integer  :: ind
        x   = log(a)

        call self%EFTLambda%precompute( x, ind, mu )

        ! compute the EFT functions:
        eft_cache%EFTOmegaV    = 0._dl
        eft_cache%EFTOmegaP    = 0._dl
        eft_cache%EFTOmegaPP   = 0._dl
        eft_cache%EFTOmegaPPP  = 0._dl
        eft_cache%EFTc         = self%EFTc%value( x, index=ind, coeff=mu )
        eft_cache%EFTLambda    = self%EFTLambda%value( x, index=ind, coeff=mu )
        eft_cache%EFTcdot      = self%EFTc%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTLambdadot = self%EFTLambda%first_derivative( x, index=ind, coeff=mu )

    end subroutine EFTCAMBADEBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the second order EFT functions at a given time.
    subroutine EFTCAMBADESecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_ADE)                            :: self          !< the base class
        real(dl), intent(in)                          :: a             !< the input scale factor
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout)  :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu
        integer  :: ind, nu_i

        real(dl) :: h0_Mpc, cs

        x   = log(a)
        h0_Mpc = eft_par_cache%h0_Mpc
        cs     = self%cs2

        ! compute the EFT functions:
        eft_cache%EFTGamma1V  = 1._dl/(4._dl*h0_Mpc**2 *a**2) * 2._dl*eft_cache%EFTc*((cs**-1) -1._dl)
        eft_cache%EFTGamma1P  = 1._dl/(4._dl*h0_Mpc**2)*(2._dl*eft_cache%EFTcdot/(eft_cache%adotoa*a**3))*((cs**-1) -1._dl)
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

    end subroutine EFTCAMBADESecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of ADE.
    subroutine EFTCAMBADEInitBackground( self, params_cache, feedback_level, success, outroot )

        implicit none

        class(EFTCAMB_ADE)                            :: self           !< the base class
        type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                       , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                       , intent(out)   :: success        !< wether the background initialization succeded or not
        character(LEN=*), optional    , intent(in)    :: outroot        !< the output root for the debug files

        ! some feedback:
        if ( feedback_level>0 ) then
            write(*,'(a)') "***************************************************************"
            write(*,'(a)') ' EFTCAMB ADE background solver'
            write(*,'(a)')
        end if

        ! initialize interpolating functions:
        call self%EFTLambda%initialize ( self%background_num_points, self%x_initial, self%x_final )
        call self%EFTc%initialize      ( self%background_num_points, self%x_initial, self%x_final )
        call self%EFTgamma1%initialize ( self%background_num_points, self%x_initial, self%x_final )

        ! debug code to print all background quantities:
        if ( DebugEFTCAMB ) then
            print*, 'EFTCAMB DEBUG ( ADE): Printing background results'
            file_debug_1%unit = 33
            file_debug_2%unit = 44
            call file_debug_1%CreateFile( TRIM(outroot)//'background_ADE_solution_1.dat' )
            call file_debug_2%CreateFile( TRIM(outroot)//'background_ADE_solution_2.dat' )
            write (file_debug_1%unit,'(a)')         '# a wADE+1 wADE_p OMEGA_ADE grho_ADE adotoa Hdot rhonu_tot*a**2/(3*adotoa**2) rhonu_ac'
            write (file_debug_2%unit,'(a)')         '# a EFTLambda EFTLambdap EFTc '
            call self%solve_background_equations ( params_cache, success=success )
            call file_debug_1%close()
            call file_debug_2%close()
        else
            call self%solve_background_equations ( params_cache, success=success )
        end if

    end subroutine EFTCAMBADEInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that solves the background equations.
    subroutine EFTCAMBADESolveBackgroundEquations( self, params_cache, success )

        implicit none

        class(EFTCAMB_ADE)                            :: self          !< the base class.
        type(TEFTCAMB_parameter_cache), intent(in)    :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters
        logical , intent(out)                         :: success       !< whether the calculation ended correctly or not


        integer  ::  i, nu_i
        real(dl) ::  wf, ac, grho_ac, grho_ADE, grho_ADEp, wADE, wADE_p
        real(dl) ::  grhom, grhorad, grhov, wADEtot,wADEtot_first , omegac, k, cs, cs_p
        real(dl) ::  a, adotoa, Hdot,Hdotdot, adotoaLCDM, h0_Mpc, grhormass_t, rhonu_tot, presnu_tot, rhonu, presnu, adotoaNoNu, grhormass_ac, rhonu_ac, presnu_ac, presnudot, presnudotdot, presnudot_tot, presnudotdot_tot

        integer  :: limit, iord(100), last , key, ier, neval
        real(dl) :: rtol, atol, t1, t2, result, alist(100), blist(100), rlist(100), elist(100), abserr

        ! 0) Initialize:
        success = .True.

        ! 1) Cosmological parameters:

        ! Derived ADE parameters:

        ac = 10._dl**(self%Log_ac)

        grhorad = params_cache%grhornomass +  params_cache%grhog !\f$ \rho_{rad} /m_0^2 \f$ total radiation background density today
        grhom =  params_cache%grhob +  params_cache%grhoc !\f$ \rho_{rad} /m_0^2 \f$ total radiation background density today

        rhonu_tot  = 0._dl
        presnu_tot = 0._dl
        presnu_ac = 0._dl
        rhonu_ac  = 0._dl
        if ( params_cache%Num_Nu_Massive /= 0) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
                rhonu  = 0._dl
                presnu = 0._dl
                grhormass_t= params_cache%grhormass(nu_i)
                grhormass_ac= params_cache%grhormass(nu_i)/ac**2
                call ThermalNuBack%rho_P(params_cache%nu_masses(nu_i),rhonu,presnu)
                rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                presnu_tot = presnu_tot + grhormass_t*presnu
                call ThermalNuBack%rho_P(ac*params_cache%nu_masses(nu_i),rhonu,presnu)
                rhonu_ac  = rhonu_ac + grhormass_ac*rhonu
                presnu_ac = presnu_ac + grhormass_ac*presnu
            end do
        end if


        grho_ac= (self%f_ac/(1._dl-self%f_ac))* (grhorad/ac**4 + grhom/ac**3 + params_cache%grhov+ rhonu_ac/ac**2)

        grhov = 3._dl*params_cache%h0_Mpc**2 - 2._dl*grho_ac/(1._dl+(1._dl/ac)**(3._dl+3._dl*self%wf)) -grhorad- grhom-rhonu_tot  ! ensures flatness (very small correction)

        h0_Mpc= params_cache%h0_Mpc


        ! 2) Initialize dqage:
        ! set-up the relative and absolute tollerances:
        rtol = 1.d-12
        atol = 1.d-16
        key = 2
        limit=100

        t1  = self%EFTLambda%x(1)
        call output( 1, t1)

        do i=1, self%background_num_points-1
            ! set the time step:
            t1 = self%EFTLambda%x(i)
            t2 = self%EFTLambda%x(i+1)
            call output( i+1, t2)
        end do

    contains
        !---------------------------------------------------------------------------------------------
        function integrand(aa)
        !Return the function to integrate for the computation of rho_ADE
            implicit none

            real(dl), intent(in)                      :: aa             !< the input scale factor.
            real(dl) integrand

            integrand = ((self%wf + 1._dl)/(1._dl+(10._dl**(self%Log_ac)/aa)**((3._dl+3._dl*self%wf)/self%p))**self%p)/aa

        end function integrand

        !---------------------------------------------------------------------------------------------
        !> Subroutine that computes the derivative of K w.r.t. x at the time x
        subroutine derivs( x )

            implicit none

            real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed

            ! 1) convert x in a:
            a  = Exp(x)

            wADE= (self%wf + 1._dl)/(1._dl+(ac/a)**((3._dl+3._dl*self%wf)/self%p))**self%p -1._dl

            wADE_p= (1._dl/a)*3._dl*(1._dl + self%wf)**2 * (ac/a)**(3._dl*(1._dl + self%wf)/self%p)/ (1._dl+ (ac/a)**(3._dl*(1._dl + self%wf)/self%p))**(1+self%p)

            ! Solve continuity equation for ADE:

            call dqage( integrand,ac,a,atol, rtol, key, limit, result, abserr, neval, ier, alist, blist, rlist, elist, iord, last )

            grho_ADE = exp(-3._dl*result)*grho_ac

            ! Compute neutrino density and adotoa

            rhonu_tot  = 0._dl
            presnu_tot = 0._dl
            if ( params_cache%Num_Nu_Massive /= 0) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates
                    rhonu  = 0._dl
                    presnu = 0._dl
                    grhormass_t= params_cache%grhormass(nu_i)/a**2
                    call ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i),rhonu,presnu)
                    rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                    presnu_tot = presnu_tot + grhormass_t*presnu
                end do
            end if

            adotoaNoNu = sqrt(grhorad/(3._dl*a**2) + grhom/(3._dl*a) + grhov*a**2/3._dl + grho_ADE*a**2/3._dl)
            adotoa= sqrt(grhorad/(3._dl*a**2) + grhom/(3._dl*a) + grhov*a**2/3._dl + grho_ADE*a**2/3._dl +rhonu_tot/3._dl)
            adotoaLCDM= sqrt(grhorad/(3._dl*a**2) + grhom/(3._dl*a) + grhov*a**2/3._dl +rhonu_tot/3._dl)

            Hdot=   -adotoa**2/2._dl - (grhorad/(3._dl*a**2)-grhov*a**2+ wADE* grho_ADE*a**2+ presnu_tot)/2._dl

            rhonu_tot        = 0._dl
            presnu_tot       = 0._dl
            presnudot_tot    = 0._dl
            presnudotdot_tot = 0._dl
            if ( params_cache%Num_Nu_Massive /= 0) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates
                    rhonu        = 0._dl
                    presnu       = 0._dl
                    presnudot    = 0._dl
                    presnudotdot = 0._dl
                    grhormass_t= params_cache%grhormass(nu_i)/a**2
                    call ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i),rhonu,presnu)
                    presnudot = ThermalNuBack%pidot(a*params_cache%nu_masses(nu_i),adotoa,presnu)
                    presnudotdot = ThermalNuBack%pidotdot(a*params_cache%nu_masses(nu_i),adotoa,Hdot,presnu,presnudot)

                    rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                    presnu_tot = presnu_tot + grhormass_t*presnu
                    rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                    presnu_tot = presnu_tot + grhormass_t*presnu
                    presnudot_tot  = presnudot_tot + grhormass_t*(presnudot -4._dl*adotoa*presnu)
                    presnudotdot_tot = presnudotdot_tot + grhormass_t*(presnudotdot &
                        & -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))
                end do
            end if

            Hdotdot= (grhom/a)*adotoa*(1._dl/6._dl) + (grhorad/a**2)*adotoa*(2._dl/3._dl)+(grho_ADE*a**2)*adotoa*(1._dl/6._dl+ wADE + 1.5_dl*wADE**2- 0.5_dl*a*wADE_p)+(grhov*a**2)*adotoa*2._dl/3._dl + adotoa*rhonu_tot/6._dl -adotoa*presnu_tot/2._dl - 0.5_dl* presnudot_tot

            cs= self%cs2
            cs_p= 0


        end subroutine

        !---------------------------------------------------------------------------------------------
        !> Subroutine that computes the Jacobian of the system. Now a dummy function.
        !! Implementing it might increase performances.
        subroutine jacobian( num_eq, x, y, ml, mu, pd, nrowpd )

            implicit none

            integer                            :: num_eq !< number of components of the Jacobian
            integer                            :: ml     !< ignored
            integer                            :: mu     !< ignored
            integer                            :: nrowpd !< ignored
            real(dl)                           :: x      !< time at which the Jacobian is computed
            real(dl), dimension(num_eq)        :: y      !< input status of the system
            real(dl), dimension(nrowpd,num_eq) :: pd     !< output Jacobian

        end subroutine jacobian

        !---------------------------------------------------------------------------------------------
        !> Subroutine that takes the solution of the JBD equations and computes additional quantities
        !! and the values of the EFT functions.
        subroutine output( ind, x)

            implicit none

            integer , intent(in)                     :: ind    !< index of the EFT functions interpolation tables to fill.
            real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed.

            logical :: is_open

            call derivs(x)
            ! 1) convert x in a:
            a  = Exp(x)

            ! 2) compute the EFT functions:
            ! NB: The Planck mass is included in grhom

            self%EFTLambda%y(ind)  =  wADE*grho_ADE*a**2 -grhov*a**2
            self%EFTLambda%yp(ind) = a**2*grho_ADE*adotoa*(a*wADE_p-3._dl*wADE*(1._dl+wADE))

            self%EFTc%y(ind)  = 0.5_dl*a**2 * (grho_ADE)*(1._dl+ wADE)
            self%EFTc%yp(ind) =  (grho_ADE)*adotoa*(-3._dl*(wADE+1._dl)**2+ a*wADE_p)*a**2/2._dl

            self%EFTgamma1%y(ind)  = 1._dl/(4._dl*h0_Mpc**2 *a**2) * 2._dl*self%EFTc%y(ind) *((cs**-1) -1._dl)
            self%EFTgamma1%yp(ind) =  1._dl/(4._dl*h0_Mpc**2)*2._dl*self%EFTc%yp(ind)/(adotoa*a**3)*((cs**-1) -1._dl)

            ! write debug quantities:
            if ( DebugEFTCAMB ) then
                inquire( unit=file_debug_1%unit, opened=is_open )
                if ( is_open ) then
                    write (file_debug_1%unit,'(20E15.5)') a,wADE+1, wADE_p, grho_ADE*a**2/(3._dl*adotoa**2),grho_ADE, adotoa, Hdot, rhonu_tot*a**2/(3*adotoa**2),rhonu_ac
                end if
                inquire( unit=file_debug_2%unit, opened=is_open )
                if ( is_open ) then
                    write (file_debug_2%unit,'(20E15.5)') a, self%EFTLambda%y(ind) , self%EFTLambda%yp(ind),self%EFTc%y(ind)
                end if
            end if
        end subroutine

    end subroutine EFTCAMBADESolveBackgroundEquations

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBADEAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_ADE)                           :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBADEAdditionalModelStability           !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBADEAdditionalModelStability = .True.


        if ( self%Log_ac > 0.0_dl ) then
            print*, 'Chose an appropriate value for parameter Log(a_c)'
            EFTCAMBADEAdditionalModelStability = .False.
            return
        end if

        if ( self%wf< 0.0_dl) then
            print*, 'Chose an appropriate value for parameter wf'
            EFTCAMBADEAdditionalModelStability = .False.
            return
        end if

        if (self%f_ac >= 1.0_dl .or. self%f_ac <= 0.0_dl) then
            print*, 'Chose an appropriate value for parameter f_ac'
            EFTCAMBADEAdditionalModelStability = .False.
            return
        end if

    end function EFTCAMBADEAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_ADE_Mod

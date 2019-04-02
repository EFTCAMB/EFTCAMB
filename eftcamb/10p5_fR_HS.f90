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

!> @file 10p5_fR_HS.f90
!! This file contains the definition of the Hu-Sawicki f(R) model.
!! Please refer to the numerical notes for details.

!----------------------------------------------------------------------------------------
!> This module contains the definition of the Hu-Sawicki f(R) model.
!! Please refer to the numerical notes for details.

!> @author Matteo Rizzato, Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_full_fR_HS

    use precision
    use IniFile
    use AMLutils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_model_full

    implicit none

    private

    public EFTCAMB_fR_HS

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of the Hu-Sawicki f(R) model.
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_fR_HS

        ! the model parameters:
        real(dl)  :: fm_hs_n        !< Hu-Sawicki f(R) model parameter n
        real(dl)  :: fm_hs_fro      !< Hu-Sawicki f(R) model parameter fR_0

        real(dl), parameter :: solver_switch = -5._dl

        ! the interpolated EFT functions that come out of the background sover:
        type(equispaced_linear_interpolate_function_1D) :: EFTOmega       !< The interpolated function Omega (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTLambda      !< The interpolated function Lambda (and derivatives).

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBHSfRReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBHSfRAllocateModelSelection      !< subroutine that allocates the model selection. For Horava this is a dummy procedure.
        procedure :: init_model_parameters           => EFTCAMBHSfRInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBHSfRInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBHSfRInitBackground               !< subroutine that initializes the background of Hu-Sawicki f(R).
        procedure :: solve_background_equations      => EFTCAMBHSfRSolveBackgroundEquations       !< subroutine that solves the Hu-Sawicki f(R) background equations.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBHSfRComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBHSfRFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBHSfRParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBHSfRParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBHSfRParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBHSfRBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBHSfRSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_adotoa                    => EFTCAMBHSfRComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBHSfRComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.

        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBHSfRAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_fR_HS

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBHSfRReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_fR_HS)       :: self   !< the base class
        type(TIniFile)             :: Ini    !< Input ini file

    end subroutine EFTCAMBHSfRReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBHSfRAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_fR_HS)       :: self !< the base class

    end subroutine EFTCAMBHSfRAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBHSfRInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_fR_HS)                          :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

        self%fm_hs_n   = array(1)
        self%fm_hs_fro = array(1)

    end subroutine EFTCAMBHSfRInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBHSfRInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_fR_HS)           :: self   !< the base class
        type(TIniFile)                 :: Ini    !< Input ini file

        self%fm_hs_n   = Ini_Read_Double_File( Ini, 'HS_n'   , 1._dl      )
        self%fm_hs_fro = Ini_Read_Double_File( Ini, 'HS_fR_0', -0.0001_dl )

    end subroutine EFTCAMBHSfRInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of Hu-Sawicki f(R).
    subroutine EFTCAMBHSfRInitBackground( self, params_cache, feedback_level, success )

        implicit none

        class(EFTCAMB_fR_HS)                         :: self           !< the base class
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                      , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                      , intent(out)   :: success        !< wether the background initialization succeded or not

        real(dl) :: fm_hs_Omegadeh2,fm_hs_Omegagammah2, fm_hs_Omeganuh2, fm_hs_Omegacdmh2
        real(dl) :: fm_hs_c2, fm_hs_c1, fm_hs_m2, fm_hs_mgamma2, fm_hs_rnge, fm_hs_step

        ! fm_hs_llim = fm_xInitial
        ! fm_hs_ulim = fm_xFinal
        ! fm_hs_nnosteps = fm_nstep+1
        ! fm_hs_z_ast_threshold = 5.d-5

        !the following is useful for debugging
        fm_hs_Omegadeh2 = params_cache%omegav*(params_cache%h0)**2.
        fm_hs_Omegagammah2 = params_cache%omegag*(params_cache%h0)**2.
        fm_hs_Omeganuh2 = params_cache%omegan*(params_cache%h0)**2.
        fm_hs_Omegacdmh2 = params_cache%omegab*(params_cache%h0)**2. +params_cache%omegac*(params_cache%h0)**2.

        fm_hs_c2 = -fm_hs_n*(6.*fm_hs_Omegadeh2/fm_hs_Omegacdmh2)*( 3.+12.*fm_hs_Omegadeh2/fm_hs_Omegacdmh2)**(-fm_hs_n-1.)*(1./fm_hs_fro)
        fm_hs_c1 = 6.*(fm_hs_Omegadeh2/fm_hs_Omegacdmh2)*fm_hs_c2
        fm_hs_m2 = ((8320.5)**(-2.))*(fm_hs_Omegacdmh2/0.13)
        fm_hs_mgamma2 = fm_hs_Omegagammah2*(1./9.)*10**(-6.)
        fm_hs_rnge = ABS(fm_hs_llim-fm_hs_ulim)
        fm_hs_step = fm_hs_rnge/real(fm_nstep)

        ! solve the background equations and store the solution:
        call self%solve_background_equations( params_cache, success=success )

        return

    end subroutine EFTCAMBHSfRInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that solves the f(R) background equations.
    subroutine EFTCAMBHSfRSolveBackgroundEquations( self, params_cache, success )

        implicit none

        class(EFTCAMB_fR_HS)                         :: self          !< the base class.
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
        logical , intent(out)                        :: success       !< whether the calculation ended correctly or not

        ! Initialization of LSODA:
        neq  = 2
        t0   = fm_hs_llim
        tfin = fm_hs_ulim
        ! set-up the relative and absolute tollerances:
        itol = 1
        rtol = 1.d-10
        atol = 1.d-14
        ! initialize task to do:
        itask  = 1
        istate = 1
        iopt   = 1
        ! initialize the work space:
        LRN = 20 + 16*NEQ
        LRS = 22 + 9*NEQ + NEQ**2
        LRW = max(LRN,LRS)
        LIS = 20 + NEQ
        LIN = 20
        LIW = max(LIS,LIN)
        ! allocate the arrays:
        allocate(rwork(LRW))
        allocate(iwork(LIW))
        allocate(y(neq))
        ! optional lsoda input:
        RWORK(5) = 0._dl  ! the step size to be attempted on the first step. The default value is determined by the solver.
        RWORK(6) = 0._dl  ! the maximum absolute step size allowed. The default value is infinite.
        RWORK(7) = 0._dl  ! the minimum absolute step size allowed. The default value is 0.
        IWORK(5) = 0  !  flag to generate extra printing at method switches. IXPR = 0 means no extra printing (the default). IXPR = 1 means print data on each switch.
        IWORK(6) = 100  !  maximum number of (internally defined) steps allowed during one call to the solver. The default value is 500.
        IWORK(7) = 0  !  maximum number of messages printed (per problem) warning that T + H = T on a step (H = step size). This must be positive to result in a non-default value.  The default value is 10.
        IWORK(8) = 0  !  the maximum order to be allowed for the nonstiff (Adams) method.  the default value is 12.
        IWORK(9) = 0  !  the maximum order to be allowed for the stiff (BDF) method.  The default value is 5.
        ! additional lsoda stuff:
        CALL XSETF(0) ! suppress odepack printing
        ! Jacobian mode: 1=fullJacobian, 2=not provided
        JacobianMode = 2


        !-----------runge-kutta-(cdm)-------------------------------------!
        fm_hs_i = 1.

        ! set initial conditions:
        ! initial time:
        fm_hs_t = fm_hs_llim
        ! guess and guess first derivative:
        fm_hs_yH(1) = fm_hs_sourceFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)/fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
        fm_hs_v(1) = fm_hs_firstderive(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
        ! initial conditions:
        fm_hs_yHtrue(1) = fm_hs_yH(1)*EXP(-3.*fm_hs_t)+fm_hs_c1/(fm_hs_c2*6.)
        fm_hs_vtrue(1) = -3.*EXP(-3.*fm_hs_t)*fm_hs_yH(1)+EXP(-3.*fm_hs_t)*fm_hs_v(1)
        ! derived quantities:
        fm_hs_eqstate(1) = -1-(fm_hs_vtrue(1)/(3*fm_hs_yHtrue(1)))
        fm_hs_rhoeff(1) = 3*fm_hs_m2*fm_hs_yHtrue(1)
        fm_hs_i = fm_hs_i+1.
        ! initial conditions for LSODA:
        y(1) = fm_hs_yHtrue(1)
        y(2) = fm_hs_vtrue(1)

        ! solve the ODE:
        DO fm_hs_ni = 2,fm_hs_nnosteps,1

            t1 = fm_hs_t
            t2 = fm_hs_llim+fm_hs_step*(fm_hs_i-1.)
            fm_hs_t = t2

            if ( .false. ) then

                !... evolving full non-linear eq. using rk4 method
                fm_hs_h1 = fm_hs_step*fm_hs_vtrue(fm_hs_ni-1)
                fm_hs_g1 = fm_hs_step*fm_hs_derivativeForV(fm_hs_llim+(fm_hs_i-2.)*fm_hs_step,fm_hs_yHtrue(fm_hs_ni-1),fm_hs_vtrue(fm_hs_ni-1),fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n)
                fm_hs_h2 = fm_hs_step*(fm_hs_vtrue(fm_hs_ni-1)+fm_hs_h1/2)
                fm_hs_g2 = fm_hs_step*fm_hs_derivativeForV(fm_hs_llim+(fm_hs_i-2.)*fm_hs_step+fm_hs_step/2.,fm_hs_yHtrue(fm_hs_ni-1)+fm_hs_h1/2,&
                    &fm_hs_vtrue(fm_hs_ni-1)+fm_hs_g1/2.,fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n)
                fm_hs_h3 = fm_hs_step*(fm_hs_vtrue(fm_hs_ni-1)+fm_hs_h2/2)
                fm_hs_g3 = fm_hs_step*fm_hs_derivativeForV(fm_hs_llim+fm_hs_step*(fm_hs_i-2.)+fm_hs_step/2.,fm_hs_yHtrue(fm_hs_ni-1)+fm_hs_h2/2,fm_hs_vtrue(fm_hs_ni-1)+&
                    &fm_hs_g2/2.,fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n)
                fm_hs_h4 = fm_hs_step*(fm_hs_vtrue(fm_hs_ni-1)+fm_hs_h3)
                fm_hs_g4 = fm_hs_step*fm_hs_derivativeForV(fm_hs_llim+(fm_hs_i-2.)*fm_hs_step+fm_hs_step,fm_hs_yHtrue(fm_hs_ni-1)+fm_hs_h3,fm_hs_vtrue(fm_hs_ni-1)+fm_hs_g3,&
                    &fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n)

                fm_hs_v(fm_hs_ni) = fm_hs_firstderive(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
                fm_hs_yH(fm_hs_ni) = fm_hs_sourceFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)/fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
                fm_hs_vtrue(fm_hs_ni) = fm_hs_vtrue(fm_hs_ni-1) + (fm_hs_g1+fm_hs_g2*2.+fm_hs_g3*2.+fm_hs_g4)/6.
                fm_hs_yHtrue(fm_hs_ni) = fm_hs_yHtrue(fm_hs_ni-1) + (fm_hs_h1+fm_hs_h2*2.+fm_hs_h3*2.+fm_hs_h4)/6.

            else if ( fm_hs_t < solver_switch ) then

                !... using the particular guess
                fm_hs_yH(fm_hs_ni) = fm_hs_sourceFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)/fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
                fm_hs_v(fm_hs_ni) = fm_hs_firstderive(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)

                fm_hs_yHtrue(fm_hs_ni) = fm_hs_yH(fm_hs_ni)*EXP(-3.*fm_hs_t)+fm_hs_c1/(fm_hs_c2*6.)
                fm_hs_vtrue(fm_hs_ni) = -3.*EXP(-3.*fm_hs_t)*fm_hs_yH(fm_hs_ni)+EXP(-3.*fm_hs_t)*fm_hs_v(fm_hs_ni)

                y(1) = fm_hs_yHtrue(fm_hs_ni)
                y(2) = fm_hs_vtrue(fm_hs_ni)

            else if ( fm_hs_t >= solver_switch ) then

                call DLSODA ( FM_HS_derivs, neq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, FM_HS_Jac, JacobianMode)

                ! check the quality of the output:
                if ( istate<0 ) then

                    fm_hs_yH(fm_hs_ni) = fm_hs_sourceFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)/fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
                    fm_hs_v(fm_hs_ni) = fm_hs_firstderive(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)

                    fm_hs_yHtrue(fm_hs_ni) = fm_hs_yH(fm_hs_ni)*EXP(-3.*fm_hs_t)+fm_hs_c1/(fm_hs_c2*6.)
                    fm_hs_vtrue(fm_hs_ni) = -3.*EXP(-3.*fm_hs_t)*fm_hs_yH(fm_hs_ni)+EXP(-3.*fm_hs_t)*fm_hs_v(fm_hs_ni)

                    y(1) = fm_hs_yHtrue(fm_hs_ni)
                    y(2) = fm_hs_vtrue(fm_hs_ni)

                    istate = 1

                end if

                fm_hs_yHtrue(fm_hs_ni) = y(1)
                fm_hs_vtrue(fm_hs_ni)  = y(2)

            end if

            fm_hs_eqstate(fm_hs_ni) = -1-(fm_hs_vtrue(fm_hs_ni)/(3*fm_hs_yHtrue(fm_hs_ni)))
            fm_hs_rhoeff(fm_hs_ni) = 3*fm_hs_m2*fm_hs_yHtrue(fm_hs_ni)
            fm_hs_i = fm_hs_i+1.

        END DO

        call FM_HS_output(success) !... calculate background and EFT functions !bh

        return

    contains

        ! ! ---------------------------------------------------------------------------------------------
        ! !> Subroutine that computes y' given y for the background of HS f(R)
        ! subroutine derivs( num_eq, x, y, ydot )
        subroutine FM_HS_derivs( num_eq, x, y, ydot )

            implicit none

            integer  :: num_eq
            real(dl) :: x
            real(dl), dimension(num_eq) :: y
            real(dl), dimension(num_eq) :: ydot

            real(dl) :: dummy

            dummy = fm_hs_derivativeForV( x, y(1), y(2), fm_hs_m2, fm_hs_c1, fm_hs_c2, fm_hs_n)

            ydot(1) = y(2)
            ydot(2) = dummy

        end subroutine FM_HS_derivs

        ! ! ---------------------------------------------------------------------------------------------
        ! !> Subroutine that computes the Jacobian of the system. Now a dummy function.
        ! !! Implementing it might increase performances.
        ! subroutine jacobian( num_eq, x, y, ml, mu, pd, nrowpd )
        subroutine FM_HS_Jac( num_eq, x, y, ml, mu, pd, nrowpd )

            implicit none

            integer  :: num_eq
            integer  :: ml ! ignored
            integer  :: mu ! ignored
            integer  :: nrowpd
            real(dl) :: x
            real(dl), dimension(num_eq) :: y
            real(dl), dimension(nrowpd,num_eq) :: pd

        end subroutine FM_HS_Jac

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that takes the solution of the background f(R) equations and computes the value of
        !! B and stores the values of the EFT functions.
        ! subroutine output( num_eq, ind, x, y, B )
        !---------------------------------------------------------------------------!
        !... calculate background and EFT functions
        subroutine FM_HS_output(success)

            implicit none

            logical, intent(inout) :: success
            integer :: fm_hs_ni
            real(dl) :: fm_hs_i
            real(dl) :: fm_hs_rhoeffprime
            real(dl) :: fm_hs_t,fm_hs_rho,fm_hs_rhop,fm_hs_Ricci,fm_hs_Hubble,fm_hs_Hprime

            !--------------------------------Hybride-Part---------------------------------!
            fm_hs_i=1.
            DO fm_hs_ni = 1,fm_hs_nnosteps,1
                fm_hs_t = fm_hs_llim+fm_hs_step*(fm_hs_i-1.)
                fm_hs_efold(fm_hs_ni) = fm_hs_t
                fm_hs_rhoeffprime = 3*fm_hs_m2*fm_hs_vtrue(fm_hs_ni)
                fm_hs_rho = fm_hs_rhofunction(fm_hs_t,fm_hs_mgamma2,fm_hs_Neff,fm_hs_Omeganuh2)
                fm_hs_rhop = fm_hs_rhoprime(fm_hs_t,fm_hs_mgamma2,fm_hs_Neff,fm_hs_Omeganuh2)
                fm_hs_Hubble = SQRT(fm_hs_m2*EXP(-3.*fm_hs_t) + fm_hs_mgamma2*EXP(-4.*fm_hs_t) + fm_hs_rho/3. +&
                    fm_hs_rhoeff(fm_hs_ni)/3.)
                fm_hs_Hprime = ((-3*fm_hs_m2)/EXP(3*fm_hs_t) - (4*fm_hs_mgamma2)/EXP(4*fm_hs_t) + fm_hs_rhop/3. + fm_hs_rhoeffprime/3.)/(2*fm_hs_Hubble)
                fm_hs_Ricci = 6*fm_hs_Hprime*fm_hs_Hubble + 12*(fm_hs_Hubble**2)
                fm_hs_Lambda(fm_hs_ni) = (1./2.)*(-(fm_hs_c1/fm_hs_c2)*fm_hs_m2 + (fm_hs_c1/(fm_hs_c2**2.))*((fm_hs_m2/fm_hs_Ricci)**(fm_hs_n+1.))*fm_hs_Ricci -&
                    fm_hs_Ricci*(-(fm_hs_c1/(fm_hs_c2**2.))*fm_hs_n*(fm_hs_m2/fm_hs_Ricci)**(1.+fm_hs_n)))
                fm_hs_Omega(fm_hs_ni) = -(fm_hs_c1/(fm_hs_c2**2.))*fm_hs_n*((fm_hs_m2/fm_hs_Ricci)**(1.+fm_hs_n))
                fm_hs_c(fm_hs_ni) = 0.
                fm_hs_i = fm_hs_i+1.
            END DO
            success = .true.

            return

        end subroutine FM_HS_output

        !---------------------------------------------------------------------------!
        function fm_hs_rhofunction(t,mgamma2,Neff,Omeganuh2)

            implicit none
            real(dl)::t,mgamma2,Neff,Omeganuh2
            real(dl) :: fm_hs_rhofunction

            IF (Omeganuh2==0.) then
                fm_hs_rhofunction = 0.
            else
                fm_hs_rhofunction = (0.6813*mgamma2*Neff*(1 + 5.434509745635465e8*(EXP(t)*Omeganuh2)**1.83)**0.5464480874316939)/EXP(4*t)
            end if

        end function fm_hs_rhofunction

        !---------------------------------------------------------------------------!
        function fm_hs_rhoprime(t,mgamma2,Neff,Omeganuh2)

            implicit none
            real(dl)::t,mgamma2,Neff,Omeganuh2
            real(dl) :: fm_hs_rhoprime

            if (Omeganuh2==0.) then
                fm_hs_rhoprime = 0.
            else
                fm_hs_rhoprime = (3.7025315D8*mgamma2*Neff*Omeganuh2*(EXP(t)*Omeganuh2)**83D-2)/&
                    (EXP(3.*t)*(1.+ 5.43451D8*(EXP(t)*Omeganuh2)**1.83)**453552D-6)&
                    - (2.7252*mgamma2*Neff*EXP(-4.*t))*((1. + 5.43451D8*(EXP(t)*Omeganuh2)**1.83)**(1./1.83))
            end if

        end function fm_hs_rhoprime

        !---------------------------------------------------------------------------!
        function fm_hs_derivativeForV(t,yh,v,m2,c1,c2,n)

            implicit none
            real(dl)::t,yh,v,yhp,m2,c1,c2,n
            real(dl)::firstpart,secondpart,thirdpart,fourthpart
            real(dl) :: fm_hs_derivativeForV

            firstpart = (c1*n*(1/(3/EXP(3*t) + 3*v + 12*yH))**(1 + n)*(-1/(2.*EXP(3*t)) - yH + (3*v + 12*yH)/6.))/c2**2
            secondpart = (-((c1*m2)/c2) + (c1*m2*(1/(3/EXP(3*t) + 3*v + 12*yH))**n)/c2**2)/(6.*m2)
            thirdpart = (2*c1*n*(1/(3/EXP(3*t) + 3*v + 12*yH))**(2 + n))/(c2**2*m2)
            fourthpart = (c1*(-1 + n)*n*(1/(3/EXP(3*t) + 3*v + 12*yH))**(2 + n))/(c2**2*m2)
            fm_hs_derivativeForV =  -4*v+(1./3.)*(9*EXP(-3*t) - (yH+firstpart + secondpart)/(m2*(EXP(-3*t)+yH)*(thirdpart+fourthpart)))

        end function fm_hs_derivativeForV

        !-------------------------------------------------------------------------------------------!
        function fm_hs_sourceFunction(m2,c1,c2,n,t)

            implicit none
            real(dl),intent(in)::m2,c1,c2,n,t
            real(dl) :: fm_hs_sourceFunction

            fm_hs_sourceFunction =3 + (3*c2)/((6*c2 + c1*EXP(3*t))*(1 + n)) +&
                (c1*EXP(3*t))/((6*c2 + c1*EXP(3*t))*(1 + n)) - &
                (2*c1**2*EXP(6*t))/(3.*c2*(6*c2 + c1*EXP(3*t))*(1 + n)) - &
                (3*c2)/((6*c2 + c1*EXP(3*t))*n*(1 + n)) - &
                (4*c1*EXP(3*t))/((6*c2 + c1*EXP(3*t))*n*(1 + n)) - &
                (4*c1**2*EXP(6*t))/(3*c2*(6*c2 + c1*EXP(3*t))*n*(1 + n))

        end function fm_hs_sourceFunction

        !-------------------------------------------------------------------------------------------!
        function fm_hs_massFunction(m2,c1,c2,n,t)

            implicit none
            real(dl),intent(in) :: m2,c1,c2,n,t
            real(dl) :: fm_hs_massFunction

            fm_hs_massFunction =-3 - (18*c2**2)/((6*c2 + c1*EXP(3*t))**2*(1 + n)) -&
                (6*c1*c2*EXP(3*t))/((6*c2 + c1*EXP(3*t))**2*(1 + n)) - &
                (5*c1**2*EXP(6*t))/((6*c2 + c1*EXP(3*t))**2*(1 + n)) + &
                (18*c2**2)/((6*c2 + c1*EXP(3*t))**2*n*(1 + n)) + &
                (6*c1*c2*EXP(3*t))/((6*c2 + c1*EXP(3*t))**2*n*(1 + n)) - &
                (4*c1**2*EXP(6*t))/((6*c2 + c1*EXP(3*t))**2*n*(1 + n)) + &
                (162*c2**3)/&
                ((1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*&
                (1 + n)) + (108*c2**4)/&
                (c1*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*&
                (6*c2 + c1*EXP(3*t))**2*n*(1 + n)) + &
                (72*c1*c2**2*EXP(3*t))/&
                ((1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*&
                (1 + n)) + (8*c1**2*c2*EXP(6*t))/&
                ((1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*&
                (1 + n))

        end function fm_hs_massFunction

        !-------------------------------------------------------------------------------------------!
        function fm_hs_fFunction(m2,c1,c2,n,t)

            implicit none
            real(dl),intent(in) :: m2,c1,c2,n,t
            real(dl) :: fm_hs_fFunction

            fm_hs_fFunction = (EXP(-3.*t)*(3.*c1*n*EXP(3.*t)+4.*c1*EXP(3.*t)+6.*c2))/(n*(n+1.)*(c1*EXP(3.*t)+6.*c2))
            fm_hs_fFunction = 2.*EXP(-3.*t) - fm_hs_fFunction
            fm_hs_fFunction = -EXP(-3.*t)*fm_hs_fFunction

        end function fm_hs_fFunction

        !-------------------------------------------------------------------------------------------!
        function fm_hs_est_z_ast(success)

            implicit none
            real(dl) :: fm_hs_est_z_ast
            integer fm_hs_ni
            real(dl) fm_hs_i,fm_hs_t,fm_hs_ratio
            logical, intent(inout) :: success

            if (ABS(fm_hs_fro) .gt. 1.d-2) then
                fm_hs_i = 1.
                DO fm_hs_ni = 1,fm_hs_nnosteps,1
                    fm_hs_t = fm_hs_llim+fm_hs_step*(fm_hs_i-1.)
                    fm_hs_yH(fm_hs_ni) = fm_hs_sourceFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)/fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
                    fm_hs_v(fm_hs_ni) = fm_hs_firstderive(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)
                    if ((fm_hs_yH(fm_hs_ni) .ne. 0.0) .and. (fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t) .ne. 0.0)) then
                        !fm_hs_ratio = fm_hs_fFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)*fm_hs_v(fm_hs_ni)
                        fm_hs_ratio = fm_hs_v(fm_hs_ni)
                        fm_hs_ratio = fm_hs_ratio/fm_hs_massFunction(fm_hs_m2,fm_hs_c1,fm_hs_c2,fm_hs_n,fm_hs_t)/fm_hs_yH(fm_hs_ni)
                        fm_hs_ratio = ABS(fm_hs_ratio)
                    endif
                    if ( fm_hs_ratio .gt. fm_hs_z_ast_threshold) then
                        fm_hs_est_z_ast = fm_hs_t
                        return
                    endif
                    fm_hs_i = fm_hs_i+1.
                END DO
                fm_hs_est_z_ast = fm_xFinal !bh
            else
                fm_hs_est_z_ast = fm_xFinal !bh
            endif
            success = .true.

        end function fm_hs_est_z_ast

        !-------------------------------------------------------------------------------------------!
        function fm_hs_firstderive(m2,c1,c2,n,t)

            implicit none
            real(dl),intent(in)::m2,c1,c2,n,t
            real(dl) :: fm_hs_firstderive

            fm_hs_firstderive = -(-27/EXP(3*t) - (6*c1*(3*c2 + 2*c1*EXP(3*t) - 3*c2*n + c1*EXP(3*t)*n))/(c2*(6*c2 + c1*EXP(3*t))*n*(1 + n)) +&
                (3*c1*(3*c2 + 2*c1*EXP(3*t))*(3*c2 + 2*c1*EXP(3*t) - 3*c2*n + c1*EXP(3*t)*n))/&
                (c2*(6*c2 + c1*EXP(3*t))**2*n*(1 + n)) +&
                (3*(3*c2 + 2*c1*EXP(3*t))*(3*c2 + 2*c1*EXP(3*t) - 3*c2*n + c1*EXP(3*t)*n))/&
                (c2*EXP(3*t)*(6*c2 + c1*EXP(3*t))*n*(1 + n)) -&
                ((3*c2 + 2*c1*EXP(3*t))*(6*c1*EXP(3*t) + 3*c1*EXP(3*t)*n))/(c2*EXP(3*t)*(6*c2 + c1*EXP(3*t))*n*(1 + n)))/&
                (3*(3/EXP(3*t) + (-108*c2**4 - 162*c1*c2**3*EXP(3*t) - 72*c1**2*c2**2*EXP(6*t) - 8*c1**3*c2*EXP(9*t) -&
                18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n - 6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n +&
                4*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n + 18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
                6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n + 5*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n)/&
                (c1*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*(1 + n)))) +&
                ((9/EXP(3*t) - ((3*c2 + 2*c1*EXP(3*t))*(3*c2 + 2*c1*EXP(3*t) - 3*c2*n + c1*EXP(3*t)*n))/&
                (c2*EXP(3*t)*(6*c2 + c1*EXP(3*t))*n*(1 + n)))*&
                (-9/EXP(3*t) - (9*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 - n)*&
                (-108*c2**4 - 162*c1*c2**3*EXP(3*t) - 72*c1**2*c2**2*EXP(6*t) - 8*c1**3*c2*EXP(9*t) -&
                18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n -&
                6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n + 4*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n +&
                18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
                6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n + 5*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n))&
                /(c1*EXP(9*t)*(6*c2 + c1*EXP(3*t))**2*(1 + n)) -&
                (6*(-108*c2**4 - 162*c1*c2**3*EXP(3*t) - 72*c1**2*c2**2*EXP(6*t) - 8*c1**3*c2*EXP(9*t) -&
                18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n -&
                6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n + 4*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n +&
                18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
                6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n + 5*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n))&
                /(EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**3*n*(1 + n)) -&
                (6*(-108*c2**4 - 162*c1*c2**3*EXP(3*t) - 72*c1**2*c2**2*EXP(6*t) - 8*c1**3*c2*EXP(9*t) -&
                18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n -&
                6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n + 4*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n +&
                18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
                6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n + 5*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n))&
                /(c1*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*(1 + n)) +&
                (-486*c1*c2**3*EXP(3*t) - 432*c1**2*c2**2*EXP(6*t) - 72*c1**3*c2*EXP(9*t) -&
                54*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n - 36*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n +&
                36*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n + 54*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
                36*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
                45*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n - 162*c1*c2**2*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 + n)*n -&
                54*c1**2*c2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 + n)*n +&
                36*c1**3*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 + n)*n +&
                162*c1*c2**2*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 + n)*n**2 +&
                54*c1**2*c2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 + n)*n**2 +&
                45*c1**3*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**(1 + n)*n**2)/&
                (c1*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*(1 + n))))/&
                (3*(3/EXP(3*t) + (-108*c2**4 - 162*c1*c2**3*EXP(3*t) - 72*c1**2*c2**2*EXP(6*t) - 8*c1**3*c2*EXP(9*t) -&
                18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n - 6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n +&
                4*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n + 18*c1*c2**2*EXP(3*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n +&
                6*c1**2*c2*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n + 5*c1**3*EXP(9*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*n)/&
                (c1*EXP(6*t)*(1/((2*c1)/c2 + 3/EXP(3*t)))**n*(6*c2 + c1*EXP(3*t))**2*n*(1 + n)))**2)

        end function fm_hs_firstderive

    end subroutine EFTCAMBHSFRSolveBackgroundEquations

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBHSfRComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_fR_HS)  :: self   !< the base class

        self%parameter_number = 2

    end subroutine EFTCAMBHSfRComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBHSfRFeedback( self, print_params )

        implicit none

        class(EFTCAMB_fR_HS)  :: self         !< the base class
        logical, optional              :: print_params !< optional flag that decised whether to print numerical values
                                                       !! of the parameters.

        logical                        :: print_params_temp

        ! print general model informations:
        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
        write(*,*)   '                 n    ='  , self%fm_hs_n
        write(*,*)   '                 fR_0 ='  , self%fm_hs_fro

        ! print the values of the parameters:
        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

    end subroutine EFTCAMBHSfRFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBHSfRParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_fR_HS) :: self   !< the base class
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
            name = 'n'
            return
        end if
        if ( i==2 ) then
            name = 'fR_0'
            return
        end if
        if ( i==0 ) then
            name = 'no_name'
            return
        end if

    end subroutine EFTCAMBHSfRParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBHSfRParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_fR_HS) :: self       !< the base class
        integer     , intent(in)      :: i          !< The index of the parameter
        character(*), intent(out)     :: latexname  !< the output latex name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names_latex.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        ! return the appropriate name:
        if ( i==1 ) then
            latexname = 'n'
            return
        end if
        if ( i==2 ) then
            latexname = 'f_{R}^0'
            return
        end if
        if ( i==0 ) then
            latexname = 'noname'
            return
        end if

    end subroutine EFTCAMBHSfRParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBHSfRParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_fR_HS) :: self   !< the base class
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
            value = self%fm_hs_n
            return
        end if
        if ( i==2 ) then
            value = self%fm_hs_fro
            return
        end if

    end subroutine EFTCAMBHSfRParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBHSfRBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_HS)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        ! compute the background EFT functions:
        eft_cache%EFTOmegaV    = 0._dl
        eft_cache%EFTOmegaP    = 0._dl
        eft_cache%EFTOmegaPP   = 0._dl
        eft_cache%EFTOmegaPPP  = 0._dl
        eft_cache%EFTc         = 0._dl
        eft_cache%EFTLambda    = 0._dl
        eft_cache%EFTcdot      = 0._dl
        eft_cache%EFTLambdadot = 0._dl

    end subroutine EFTCAMBHSfRBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the second order EFT functions at a given time.
    subroutine EFTCAMBHSfRSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_HS)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: psi, psiprime, psiprimeprime, a2, phip1, phip2, phip3, c3

        ! compute the second order EFT functions:
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

    end subroutine EFTCAMBHSfRSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H.
    subroutine EFTCAMBHSfRComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_HS)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp, a2, Omega_tot
        integer     :: nu_i

        ! temp =
        ! eft_cache%adotoa = sqrt( temp )

    end subroutine EFTCAMBHSfRComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBHSfRComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_HS)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp, a2, Omega_tot, Omega_tot_prime, Omega_tot_primeprime, Omega_phi0
        integer     :: nu_i

        eft_cache%Hdot    = 0._dl
        eft_cache%Hdotdot = 0._dl

    end subroutine EFTCAMBHSfRComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBHSfRAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_HS)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBHSfRAdditionalModelStability       !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBHSfRAdditionalModelStability = .True.

    end function EFTCAMBHSfRAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_full_fR_HS

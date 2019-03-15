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

!> @file 10p3_Kmouflage.f90
!! This file contains

!----------------------------------------------------------------------------------------
!> This module contains

!> @author Giampaolo Benevento, Marco Raveri

module EFTCAMB_Kmouflage_Mod

    use precision
    use IniFile
    use AMLutils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_model_full
    use equispaced_linear_interpolation_1D
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_constant_parametrization_1D

    implicit none

    private

    public EFTCAMB_Kmouflage

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of K-mouflage models. The reference papers are: arXiv:1809.09958 arXiv:1509.00611
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_Kmouflage


        ! the model selection flag:
        logical   :: Kmimic              !< Selects whether to consider the version of K-mouflage that mimic a LCDM background expansion.
                                                    !! Notice that when this is used the parameters alphaU and gammaU are ignored.
        ! the model parameters:
        real(dl)  :: alphaU      !< first Kmouflage model parameter.
        real(dl)  :: gammaU      !< second Kmouflage model parameter.
        real(dl)  :: m           !< third Kmouflage model parameter.
        real(dl)  :: eps2_0      !< fourth Kmouflage model parameter.
        real(dl)  :: gammaA      !< fifth Kmouflage model parameter.

        ! the interpolated EFT functions that come out of the background sover:
        type(equispaced_linear_interpolate_function_1D) :: EFTOmega       !< The interpolated function Omega (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTLambda      !< The interpolated function Lambda (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTc           !< The interpolated function c (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTgamma1      !< The interpolated function gamma1 (and derivatives).

        ! some designer parameters:
        integer  :: background_num_points = 1000                          !< Number of points sampled by the designer code.
        real(dl) :: x_initial             = log(10._dl**(-9._dl))         !< log(a start)
        real(dl) :: x_final               = 0.0_dl                        !< log(a final)

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBKmouflageReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBKmouflageAllocateModelSelection      !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBKmouflageInitModelParameters         !< subroutine that initializes the model
        procedure :: init_model_parameters_from_file => EFTCAMBKmouflageInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBKmouflageComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBKmouflageFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBKmouflageParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBKmouflageParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBKmouflageParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBKmouflageBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBKmouflageSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBKmouflageInitBackground             !< subroutine that initializes the background of Kmouflage.
        procedure :: solve_background_equations      => EFTCAMBKmouflageSolveBackgroundEquations   !< subroutine that solves the background equations.

        ! stability procedures:
        procedure :: additional_model_stability      => EFTCAMBKmouflageAdditionalModelStability   !< function that computes model specific stability requirements.

    end type EFTCAMB_Kmouflage
 
 ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBKmouflageReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Kmouflage)    :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        ! read model selection flags:
        self%Kmimic = Ini_Read_Logical_File( Ini, 'Kmimic', .False. )


    end subroutine EFTCAMBKmouflageReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection. 
    subroutine EFTCAMBKmouflageAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_Kmouflage) :: self !< the base class

        ! nothing to be done here.

    end subroutine EFTCAMBKmouflageAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBKmouflageInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_Kmouflage)                               :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

        real(dl), allocatable, dimension(:)                    :: temp
        integer :: num_params_function, num_params_temp, i

        self%alphaU  = array(1)
        self%gammaU  = array(2)
        self%m       = array(3)
        self%eps2_0  = array(4)
        self%gammaA  = array(5)

    end subroutine EFTCAMBKmouflageInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBKmouflageInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Kmouflage)       :: self   !< the base class
        type(TIniFile)                 :: Ini    !< Input ini file

        self%alphaU = Ini_Read_Double_File( Ini, 'alphaU', 0.2_dl   )
        self%gammaU = Ini_Read_Double_File( Ini, 'gammaU', 1._dl    )
        self%m      = Ini_Read_Double_File( Ini, 'm'     , 3._dl    )
        if ( self%Kmimic ) then
        self%eps2_0 = Ini_Read_Double_File( Ini, 'eps2_0' , 0.01_dl )
        else
        self%eps2_0 = Ini_Read_Double_File( Ini, 'eps2_0' , -0.01_dl )
        end if
        self%gammaA = Ini_Read_Double_File( Ini, 'gammaA', 0.2_dl   )

    end subroutine EFTCAMBKmouflageInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBKmouflageComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_Kmouflage)  :: self   !< the base class      

            self%parameter_number = 5

    end subroutine EFTCAMBKmouflageComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBKmouflageFeedback( self, print_params )

        implicit none

        class(EFTCAMB_Kmouflage)     :: self         !< the base class
        logical, optional            :: print_params !< optional flag that decised whether to print numerical values
                                                     !! of the parameters.

        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        if ( self%Kmimic ) then
            write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number-2
            write(*,'(a)')  '   K-mouflage with LCDM expansion history, parameters related to the function U are ignored'
        else
            write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
        end if

        ! print model informations:
        write(*,*)
        write(*,'(a24,F12.6)') '    alphaU       =', self%alphaU
        write(*,'(a24,F12.6)') '    gammaU       =', self%gammaU
        write(*,'(a24,F12.6)') '    m            =', self%m
        write(*,'(a24,F12.6)') '    eps2_0       =', self%eps2_0
        write(*,'(a24,F12.6)') '    gammaA       =', self%gammaA

    end subroutine EFTCAMBKmouflageFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBKmouflageParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_Kmouflage)    :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! other parameters:
        else if ( i==1 ) then
            name = TRIM('alphaU')
            return
        else if ( i==2 ) then
            name = TRIM('gammaU')
        else if ( i==3 ) then
            name = TRIM('m')
        else if ( i==4 ) then
            name = TRIM('eps2_0')
        else if ( i==5 ) then
            name = TRIM('gammaA')
            return
        end if

    end subroutine EFTCAMBKmouflageParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBKmouflageParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_Kmouflage)    :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! other parameters:
        else if ( i==1 ) then
            latexname = TRIM('\alpha_{U}')
            return
        else if ( i==2 ) then
            latexname = TRIM('\gamma_{U}')
            return
        else if ( i==3 ) then
            latexname = TRIM('m')
            return
        else if ( i==4 ) then
            latexname = TRIM('\epsilon_{2,0}')
            return
        else if ( i==5 ) then
            latexname = TRIM('\gamma_{A}')
            return
        end if

    end subroutine EFTCAMBKmouflageParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBKmouflageParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_Kmouflage)    :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! other parameters:
        else if ( i==1 ) then
            value = self%alphaU
            return
        else if ( i==2 ) then
            value = self%gammaU
            return
        else if ( i==3 ) then
            value = self%m
            return
        else if ( i==4 ) then
            value = self%eps2_0
            return
        else if ( i==5 ) then
            value = self%gammaA
            return
        end if

    end subroutine EFTCAMBKmouflageParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBKmouflageBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Kmouflage)                     :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu, grhom, grhorad
        integer  :: ind

        x   = log(a)

        call self%EFTOmega%precompute( x, ind, mu )

        ! compute the EFT functions:
        eft_cache%EFTOmegaV    = self%EFTOmega%value( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaP    = self%EFTOmega%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaPP   = self%EFTOmega%second_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaPPP  = self%EFTOmega%third_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTc         = self%EFTc%value( x, index=ind, coeff=mu )
        eft_cache%EFTLambda    = self%EFTLambda%value( x, index=ind, coeff=mu )
        eft_cache%EFTcdot      = self%EFTc%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTLambdadot = self%EFTLambda%first_derivative( x, index=ind, coeff=mu )

    end subroutine EFTCAMBKmouflageBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the second order EFT functions at a given time.
    subroutine EFTCAMBKmouflageSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Kmouflage)                     :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu
        integer  :: ind

        x   = log(a)
        call self%EFTgamma1%precompute( x, ind, mu )

        ! compute the EFT functions:
        eft_cache%EFTGamma1V  = self%EFTgamma1%value( x, index=ind, coeff=mu )
        eft_cache%EFTGamma1P  = self%EFTgamma1%first_derivative( x, index=ind, coeff=mu )
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

    end subroutine EFTCAMBKmouflageSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of Kmouflage.
    subroutine EFTCAMBKmouflageInitBackground( self, params_cache, feedback_level, success )

        implicit none

        class(EFTCAMB_Kmouflage)                     :: self           !< the base class
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                      , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                      , intent(out)   :: success        !< wether the background initialization succeded or not

        ! some feedback:
        if ( feedback_level>0 ) then
            write(*,'(a)') "***************************************************************"
            write(*,'(a)') ' EFTCAMB Kmouflage background solver'
            write(*,'(a)')
        end if

        ! protect against massive neutrinos that are not yet implemented:
        if ( params_cache%Num_Nu_Massive /= 0 ) then
            write(*,*) 'Massive neutrinos are not implemented for K-mouflage models'
            write(*,*) 'The value of Num_Nu_Massive should be 0'
            call MpiStop('EFTCAMB error')
        end if

        ! initialize interpolating functions:
        call self%EFTOmega%initialize  ( self%background_num_points, self%x_initial, self%x_final )
        call self%EFTLambda%initialize ( self%background_num_points, self%x_initial, self%x_final )
        call self%EFTc%initialize      ( self%background_num_points, self%x_initial, self%x_final )
        call self%EFTgamma1%initialize ( self%background_num_points, self%x_initial, self%x_final )

        ! solve the background equations and store the solution:

        if ( DebugEFTCAMB ) then
            print*, 'EFTCAMB DEBUG ( Kmouflage background solver ): Printing solution to background equations.'
            call CreateTxtFile( './debug_kmouflage_background_solution.dat', 33 )
            call CreateTxtFile( './debug_kmouflage_background_eft.dat', 44 )
        end if

        call self%solve_background_equations ( params_cache, success=success )

        if ( DebugEFTCAMB ) then
            close(33)
            close(44)
        end if

    end subroutine EFTCAMBKmouflageInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that solves the Klein Gordon equation to get the kinetic function K.
    subroutine EFTCAMBKmouflageSolveBackgroundEquations( self, params_cache, success )

        implicit none

        class(EFTCAMB_Kmouflage)                     :: self          !< the base class.
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
        logical , intent(out)                        :: success       !< whether the calculation ended correctly or not


        integer, parameter :: num_eq = 2                              !<  Number of equations

        real(dl) :: y(num_eq), ydot(num_eq)
        
        real(dl) :: Omegam_EFT, Omegavac_EFT, OmegaMassiveNu_EFT, OmegaGamma_EFT, OmegaNu_EFT, grhom, grhorad, OmegaPhi_0, Omegam_L
        real(dl) :: Omegarad_EFT, EquivalenceScale, a, adotoa, Hdot
        real(dl) :: U, Up, Upp, Uppp, Av, Ap, App, Appp, Apppp, Mnorm, sqrt_chi, chi, chip, chipp, DlogUoverDloga, eps2, eps2p,eps2pp, nu_A, alpha_A, dKdchi, d2Kdchi2,  d3Kdchi3, ufunc, vfunc, ufuncp,ufuncpp, vfuncp, vfuncpp, eps1, p_phi, rho_phi, OmegaPhi

        integer  :: itol, itask, istate, iopt, LRN, LRS, LRW, LIS, LIN, LIW, JacobianMode, i
        real(dl) :: rtol, atol, t1, t2
        real(dl), allocatable :: rwork(:)
        integer , allocatable :: iwork(:)
        real(dl) :: Unormalization, U1, K , Kp, adotoaL, HdotL, Hdot1

        ! 0) Initialize:
        success = .True.

        ! 1) Cosmological parameters:

        ! Derived K-mouflage parameters:

        nu_A    = 3._dl*(self%m-1._dl)/(2._dl*self%m-1._dl)                               !Eq (6.3) of 1509.00611

        alpha_A = -self%eps2_0*(self%gammaA+1._dl)/(nu_A*self%gammaA)

         
        grhorad = params_cache%grhornomass +  params_cache%grhog !\f$ \rho_{rad} /m_0^2 \f$ total radiation background density today
        
        Omegam_EFT         = params_cache%omegab + params_cache%omegac
        OmegaMassiveNu_EFT = params_cache%omegan
        OmegaGamma_EFT     = params_cache%omegag
        OmegaNu_EFT        = params_cache%omegar
        Omegarad_EFT       = OmegaGamma_EFT + OmegaNu_EFT
        Omegam_L= Omegam_EFT/(1._dl-self%eps2_0)+(2._dl*self%eps2_0*self%gammaA*(1._dl+ 2._dl*Omegarad_EFT-nu_A)+4._dl*self%eps2_0*(1._dl+ Omegarad_EFT))/((3._dl+3._dl*self%gammaA)*(1._dl-self%eps2_0))+self%eps2_0/(1._dl-self%eps2_0)
        
        if ( self%Kmimic ) then
            OmegaPhi_0 = (1._dl - self%eps2_0)**2 - (Omegam_EFT + Omegarad_EFT) 
            Omegavac_EFT       = 1._dl - Omegam_L - Omegarad_EFT
            grhom   = (params_cache%grhob + params_cache%grhoc)   !\f$ \sum_m\rho_m /m_0^2 \f$ total matter background density today
        else 
            Omegavac_EFT       = 1._dl - Omegam_EFT - Omegarad_EFT
            OmegaPhi_0 = Omegavac_EFT*(1._dl -self%eps2_0)**2 -(2._dl*self%eps2_0 - self%eps2_0**2)&
                & *(Omegam_EFT +Omegarad_EFT)
            grhom   = params_cache%grhob + params_cache%grhoc        !\f$ \sum_m\rho_m /m_0^2 \f$ total matter background density today
        end if 

        
        EquivalenceScale = (Omegarad_EFT+ OmegaMassiveNu_EFT)/Omegam_EFT

        ! 2) Set initial conditions:
        if ( self%Kmimic ) then
            y(1)  = (3._dl*(1._dl + self%gammaA)*(Omegam_L - Omegam_EFT)- self%eps2_0*(4._dl+3._dl*Omegam_L+4._dl*Omegarad_EFT+self%gammaA*(2._dl-2._dl*nu_A+3._dl*Omegam_L+4._dl*Omegarad_EFT)))&
                      &/(6._dl*(1._dl + self%gammaA)*(1._dl +(-2._dl+ self%eps2_0)*self%eps2_0  -Omegam_EFT-Omegarad_EFT)) !initial condition are (u[1]+v[1])/2u[1]
            y(2)= 0
        else 
            y(1)  = - 1._dl   !initial condition are imposed at a=1
            y(2)= 0
        end if 


        ! 3) Initialize DLSODA:
        ! set-up the relative and absolute tollerances:
        itol = 1
        rtol = 1.d-12
        atol = 1.d-16
        ! initialize task to do:
        itask  = 1
        istate = 1
        iopt   = 1
        ! initialize the work space:
        LRN = 20 + 16*num_eq
        LRS = 22 + 9*num_eq + num_eq**2
        LRW = max(LRN,LRS)
        LIS = 20 + num_eq
        LIN = 20
        LIW = max(LIS,LIN)
        ! allocate the arrays:
        allocate(rwork(LRW))
        allocate(iwork(LIW))
        ! optional lsoda input:
        RWORK(5) = 0._dl  ! the step size to be attempted on the first step. The default value is determined by the solver.
        RWORK(6) = 0._dl  ! the maximum absolute step size allowed. The default value is infinite.
        RWORK(7) = 0._dl  ! the minimum absolute step size allowed. The default value is 0.
        IWORK(5) = 0      ! flag to generate extra printing at method switches. IXPR = 0 means no extra printing (the default). IXPR = 1 means print data on each switch.
        IWORK(6) = 100    ! maximum number of (internally defined) steps allowed during one call to the solver. The default value is 500.
        IWORK(7) = 0      ! maximum number of messages printed (per problem) warning that T + H = T on a step (H = step size). This must be positive to result in a non-default value.  The default value is 10.
        IWORK(8) = 0      ! the maximum order to be allowed for the nonstiff (Adams) method.  the default value is 12.
        IWORK(9) = 0      ! the maximum order to be allowed for the stiff (BDF) method.  The default value is 5.
        ! additional lsoda stuff:
        CALL XSETF(0) ! suppress odepack printing
        ! Jacobian mode: 1=fullJacobian, 2=not provided
        JacobianMode = 2

        ! store the first value of the EFT functions:
        t1  = self%EFTOmega%x(self%EFTOmega%num_points)
        call output( num_eq, self%EFTOmega%num_points, t1, y )

        ! 4) Integrate dK/dx:
        do i= self%EFTOmega%num_points, 2, -1
            ! set the time step:
            t1 = self%EFTOmega%x(i)
            t2 = self%EFTOmega%x(i-1)
            call DLSODA ( derivs, num_eq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)

            if ( istate < 0 ) then
                success = .False.
                return
            end if
            call output( num_eq, i-1, t2, y )
        end do

    contains

        !---------------------------------------------------------------------------------------------
        !> Subroutine that computes the derivative of K w.r.t. x at the time x
        subroutine derivs( num_eq, x, y, ydot )

            implicit none

            integer , intent(in)                     :: num_eq !< number of equations in the ODE system
            real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed
            real(dl), intent(in), dimension(num_eq)  :: y      !< Imput status of the system
            real(dl), intent(out), dimension(num_eq) :: ydot   !< value of the derivative of K at time x

            ! 1) convert x in a:
            a  = Exp(x)
 
            ! 1) Define Kmouflage functions and compute derivatives
            Av      = 1._dl + alpha_A- alpha_A*(a*(self%gammaA+1._dl)/(self%gammaA+a))**nu_A  !Eq (6.2) of 1509.00611

            eps2= -(a/Av)*(alpha_A*self%gammaA*nu_A*(a*(self%gammaA+1._dl)/(self%gammaA+a))**(nu_A) ) &
                & / (a**2 + a*self%gammaA)    !Eq (2.7) of 1509.00611

            Ap= eps2*Av/a   !First derivative of $\bar{A}(a)$ w.r.t. "a"
        
            App= Av*eps2*(-2._dl*a+self%gammaA*(-1._dl+nu_A))/((a+self%gammaA)*a**2) !Second derivative of $\bar{A}(a)$ w.r.t. "a"
            
            Appp= -(1._dl/(a**3 *(a + self%gammaA)**3))* alpha_A* self%gammaA*(a*(1._dl +self%gammaA) &
                & / (a + self%gammaA))**nu_A *(6._dl* a**2 - 6._dl *a*self%gammaA*(-1._dl+nu_A) + self%gammaA**2 &
                & * (-2._dl + nu_A)* (-1._dl + nu_A)) * nu_A    !Third derivative of $\bar{A}(a)$ w.r.t. "a"

            Apppp= (self%eps2_0*((a*(1._dl + self%gammaA))/(a + self%gammaA))**(1._dl + nu_A)*(-24._dl*a**3 + 36._dl*a**2* self%gammaA*(-1._dl + nu_A) &
                 & - 12._dl*a* self%gammaA**2*(2._dl - 3._dl*nu_A + nu_A**2) + self%gammaA**3*(-6._dl + 11._dl*nu_A - 6._dl*nu_A**2 + nu_A**3)))/(a**5*(a + self%gammaA)**3)
    
            eps2p= (-(a*Ap**2) + Av*(Ap + a*App))/Av**2

            eps2pp= (2._dl*a*Ap**3 - Av*Ap*(2._dl*Ap+ 3._dl*a*App) + Av**2*(2._dl*App + a*Appp))/Av**3


            if ( self%Kmimic ) then
                Mnorm= OmegaPhi_0/Omegam_EFT
            else 
                Mnorm= OmegaPhi_0/Omegam_EFT + self%eps2_0/(-3._dl*self%eps2_0 + (2._dl*self%alphaU + &
                    &(3._dl+4._dl* sqrt(EquivalenceScale))*(1._dl+self%gammaU)*(log(1._dl+self%gammaU))**2)/(2._dl*(1._dl &
                    & + self%gammaU)*log(1._dl+self%gammaU)*(self%alphaU+(1._dl+sqrt(EquivalenceScale))&
                    &*log(1._dl+self%gammaU))))    !Eq (5.14) of 1509.00611
            end if 
            
            if ( .not. self%Kmimic ) then
                ! Part of the normalization factor of Eq (6.1) of 1509.00611:
                U1=(log(self%gammaU+1))/((sqrt(EquivalenceScale)+1)*log(self%gammaU+1)+self%alphaU) 
 
                Unormalization= U1*sqrt(-Mnorm* 2._dl*(-3._dl*self%eps2_0+ (2._dl*self%alphaU  &
                    & + (3._dl+4._dl* sqrt(EquivalenceScale))*(1._dl+self%gammaU)*(log(1._dl+self%gammaU))**2)/(2._dl*(1._dl &
                    & + self%gammaU)*log(1._dl+self%gammaU)*(self%alphaU+(1._dl+sqrt(EquivalenceScale))&
                    & * log(1._dl+self%gammaU))))/(self%eps2_0)) !  Unormalization= U1/sqrt_chi(a=1) -> Kp(a=1)=1 
            
                U=(a**2*log(self%gammaU+a)/((sqrt(EquivalenceScale)+sqrt(a))*log(self%gammaU+a)+self%alphaU*a**2)) &
                    & / Unormalization                    !Eq (6.1) of 1509.00611 normalized as K'(a=1)= 1

                DlogUoverDloga= (2._dl*self%alphaU*a**3+(3._dl*sqrt(a)+4._dl* sqrt(EquivalenceScale))*(a+self%gammaU)*&
                    &(log(a+self%gammaU))**2)/(2._dl*(a+self%gammaU)*log(a+self%gammaU)*&
                    &(self%alphaU*a**2+(sqrt(a)+sqrt(EquivalenceScale))*log(a+self%gammaU)))
       
                sqrt_chi= -(eps2*Av**4)/(2._dl*U*(-3._dl*eps2+DlogUoverDloga)* Mnorm)    !Eq (5.8) of 1509.00611
              
                chi= sqrt_chi**2
        
                Up=  DlogUoverDloga*U/a !First derivative of "U(a)" w.r.t. "a"
        
                Upp= (-4._dl*(a**4)* self%alphaU* (2._dl* sqrt(a) + 2._dl*sqrt(EquivalenceScale)  + a**2 * self%alphaU) &
                    & + log(a + self%gammaU) * (4._dl* (a**3)* self%alphaU* (5._dl* a**(1.5) + 7._dl* a* sqrt(EquivalenceScale) &
                    & + 6._dl *sqrt(a) * self%gammaU + 8._dl* sqrt(EquivalenceScale)* self%gammaU)+ ((a + self%gammaU)**2)  &
                    & * log(a + self%gammaU)* (- 3._dl* (a**2)* (5._dl* sqrt(a) + 8._dl* sqrt(EquivalenceScale))* self%alphaU &
                    & + (3._dl* a + 9._dl* sqrt(a) * sqrt(EquivalenceScale) +8._dl*EquivalenceScale) * log(a + self%gammaU)))) &
                    & /(4._dl*((a +self%gammaU)**2)*(a**2 * self%alphaU+(sqrt(a)+ sqrt(EquivalenceScale))* log(a + self%gammaU))**3 *Unormalization)
                    !Second derivative of "U(a)" w.r.t. "a"

                Uppp= (8._dl* a**(4.5_dl)* self%alphaU* (12._dl *sqrt(a)* sqrt(EquivalenceScale) + 6._dl* EquivalenceScale &
                    & + 15._dl* a**(2.5_dl)* self%alphaU + 18._dl * a**2 * sqrt(EquivalenceScale)* self%alphaU &
                    & + 2._dl* a**4 * self%alphaU**2 + 9._dl* a**(1.5_dl)* self%alphaU * self%gammaU   &
                    & + 6._dl* a* (1._dl + 2._dl* sqrt(EquivalenceScale) * self%alphaU * self%gammaU)) &
                    & + log(a + self%gammaU) * (-4._dl* a**(3.5_dl)* self%alphaU*(55._dl*a**(3.5_dl)*self%alphaU &
                    & + 88._dl* a**3 * sqrt(EquivalenceScale)* self%alphaU  + 84._dl* sqrt(a)* sqrt(EquivalenceScale) &
                    & * self%gammaU + 48._dl *EquivalenceScale *self%gammaU + 108._dl* a**(2.5_dl)*  self%alphaU &
                    & * self%gammaU + 24._dl* a**2 *(1._dl + 7._dl*sqrt(EquivalenceScale)*self%alphaU*self%gammaU)&
                    & + 15._dl* a**(1.5_dl)* (4._dl *sqrt(EquivalenceScale) + 3._dl* self%alphaU* self%gammaU**2)&
                    & + 36._dl* a* (EquivalenceScale + self%gammaU + 2._dl*sqrt(EquivalenceScale) *self%alphaU &
                    & * self%gammaU**2)) + log(a + self%gammaU)* (a**(2.5_dl)* self%alphaU* (a**2 *(88._dl* a &
                    & + 260._dl* sqrt(a) * sqrt(EquivalenceScale)+ 208._dl* EquivalenceScale + 105._dl* a**(2.5_dl)* self%alphaU &
                    & + 192._dl* a**2 *sqrt(EquivalenceScale)*self%alphaU)+ 3._dl*a*(72._dl* a+ 208._dl*sqrt(a)*sqrt(EquivalenceScale) &
                    & + 160._dl* EquivalenceScale + 105._dl* a**(2.5_dl) *self%alphaU + 192._dl* a**2 *sqrt(EquivalenceScale) &
                    & * self%alphaU)* self%gammaU + 9._dl*(16._dl*a+ 44._dl* sqrt(a)* sqrt(EquivalenceScale) + 32._dl &
                    & * EquivalenceScale+35._dl* a**(2.5_dl)*self%alphaU+ 64._dl*a**2 *sqrt(EquivalenceScale)*self%alphaU) &
                    & * self%gammaU**2 + 3._dl* a*(35._dl* sqrt(a) + 64._dl* sqrt(EquivalenceScale))* self%alphaU &
                    & * self%gammaU**3)+ 3._dl* (a + self%gammaU)**3 *log(a + self%gammaU)*(-2._dl &
                    & * a**(1.5_dl)* (10._dl* a + 33._dl* sqrt(a)*sqrt(EquivalenceScale)+ 32._dl* EquivalenceScale)* self%alphaU &
                    & - (a + 4._dl*sqrt(a)*sqrt(EquivalenceScale)+ 5._dl* EquivalenceScale) *log(a + self%gammaU)))))/(8._dl &
                    & * sqrt(a)*(a + self%gammaU)**3 *(a**2 * self%alphaU + (sqrt(a) + sqrt(EquivalenceScale)) &
                    & * log(a + self%gammaU))**4 * Unormalization) !Third derivative of "U(a)" w.r.t. "a"

                chip=  Av**7 * Ap* (12._dl* U* Ap**3 + Av*(-6._dl* Ap**2 *Up -Av* Up* App + Av* Ap *Upp)) &
                    & /(2._dl*  Mnorm**2 * (3._dl* U*Ap - Av* Up)**3)
            
                chipp= (1._dl/(2._dl * Mnorm**2 * (-3._dl* U* Ap + Av *Up)**4)) * Av**6 * (36._dl * U**2 *(7._dl* Ap**6 + Av* Ap**4 *App)+ 3._dl &
                    & * Av * U * Ap * (-88._dl* Ap**4 *Up + 2._dl* Av**2 *Up * App**2 + 15._dl* Av* Ap**3 *Upp - Av**2 *Ap* (2._dl* App* Upp &
                    & + Up *Appp)+Av*Ap**2 *(-25._dl * Up* App + Av*Uppp))+Av**2 *(72._dl* Ap**4 *Up**2 + Av**2 *Up**2 *App**2-27._dl* Av *Ap**3 &
                    & * Up* Upp + Av**2 *Ap* Up* (-4._dl *App* Upp+Up* Appp)+ Av* Ap**2 *(33._dl *Up**2 *App+ 3._dl * Av *Upp**2- Av *Up* Uppp)))
           
                dKdchi= U/(sqrt_chi * a**3)
             
                d2Kdchi2= - ((Av**4) * Ap / Mnorm + 6._dl* a**2 * Av* chi* dKdchi - 6._dl* a**3 * chi* Ap* dKdchi + a**3 * Av &
                    & * chip* dKdchi)/(2._dl * a**3 * Av* chi* chip)

                d3Kdchi3= (-3._dl* a* Av**4 * Ap**2 *chi *chip / Mnorm - 6._dl* a**4 *chi**2 *Ap**2 * chip * dKdchi &
                    & + Av**5 *(chip* (Ap* (3._dl* chi + a* chip) -a* chi* App) +  a* chi* Ap* chipp)/Mnorm + 6._dl &
                    & * a**4 * Av* chi**2 *(dKdchi *(chip * App- Ap* chipp) + Ap*chip**2 * d2Kdchi2)+ a**2 * Av**2 * (dKdchi &
                    & * (a**2 * chip**3 + 6._dl *chi**2 *(chip +a*chipp))-a*chi* chip**2 *(6._dl *chi+a*chip) &
                    & * d2Kdchi2))/(2._dl* a**4 * Av**2 *chi**2 *chip**3) !OK
   

                Kp=dKdchi*chip

                ! 2) Get the derivatives of K w.r.t. x:
                ydot(1)= dKdchi * chip*a    !Eq (5.7) of 1509.00611

                K=y(1) 

                ! 3) Compute Hubble parameter and its derivative w.r.t conformal time:

                adotoa= a *params_cache%h0_Mpc *(Av/(1-eps2))*sqrt(Omegam_EFT/a**3 + Omegarad_EFT/a**4 + OmegaPhi_0* Av**(-4)*(K-2._dl*chi*dKdchi) &
                    & /(-1._dl-2._dl* (U1/Unormalization)**2))

                ydot(2) = - sqrt( 2._dl*chi* Mnorm*grhom)* a/(Av*adotoa) 

                Hdot= (params_cache%h0_Mpc)**2*(2._dl*a*(-1._dl-2._dl* (U1/Unormalization)**2)*(a*Omegam_EFT + Omegarad_EFT)*Av**4 *(-1._dl + eps2)*Ap + 2._dl*a**5* OmegaPhi_0*(-1._dl + eps2)*(2._dl*chi*dKdchi - K) &                   
                    & * Ap -(-1._dl-2._dl* (U1/Unormalization)**2)* Av**5 * ((a*Omegam_EFT + 2._dl *Omegarad_EFT)*(-1._dl + eps2) + 2._dl*a*(a*Omegam_EFT + Omegarad_EFT)*eps2p) & 
                    & + a**4 * OmegaPhi_0 *Av*(2._dl*K*(-1._dl + eps2-a*eps2p)+chi*(-2._dl*a*(-1._dl+eps2)*d2Kdchi2*chip + 4._dl*dKdchi*(+1._dl - eps2+a*eps2p))+a*(1._dl-eps2)*(2._dl*dKdchi*chip -Kp )))/(2._dl*a**2 *(-1._dl-2._dl* (U1/Unormalization)**2) *Av**3 *(-1._dl + eps2)**3)

                adotoaL = a * params_cache%h0_Mpc * sqrt(Omegam_EFT/a**3 + Omegarad_EFT/a**4 + Omegavac_EFT)
       
                HdotL = -1.5_dl*a**2 * params_cache%h0_Mpc**2 * (Omegam_EFT/a**3 + 4._dl*Omegarad_EFT/(3._dl*a**4)) + adotoaL**2

            else
            
            !Equations (2.15-2.16) of 1809.09958 and their derivatives w.r.t. to the scale factor are rewritten in terms of ufunc and vfunc where vfunc= A^4 *p_phi and ufunc= A^4 *rho_phi
 
                ufunc= Av**2 *(1._dl-eps2)**2 * (Omegam_L/a**3 + Omegarad_EFT/a**4 + Omegavac_EFT) -Av**4 * (Omegam_EFT/a**3 + Omegarad_EFT/a**4)

                ufuncp= (-(3._dl*a*Omegam_L + 4._dl*Omegarad_EFT)*Av**2 +(3._dl*a*Omegam_EFT + 4._dl*Omegarad_EFT)*Av**4 - 4._dl*a*(a*Omegam_EFT + Omegarad_EFT)*Av**3 *Ap &
                      & - 2._dl*a*Av*((-3._dl*a*Omegam_L - 4._dl*Omegarad_EFT)*Ap + a*(a**4 *Omegavac_EFT + a*Omegam_L + Omegarad_EFT)*App) & 
                      & + a**2* Ap*((-3._dl*a*Omegam_L - 4._dl*Omegarad_EFT)*Ap  + 2._dl*a*(a**4 *Omegavac_EFT + a*Omegam_L + Omegarad_EFT)*App))/a**5

                ufuncpp= (2._dl*(-2._dl*(3._dl*a*Omegam_EFT + 5._dl*Omegarad_EFT)*Av**4 + 2._dl*Av**2*(3._dl*a*Omegam_L + 5._dl*Omegarad_EFT &
                      & - 3._dl*a**2*(a*Omegam_EFT + Omegarad_EFT)*Ap**2) + 2._dl*a*Av**3*((6._dl*a*Omegam_EFT + 8._dl*Omegarad_EFT)*Ap - a*(a*Omegam_EFT + Omegarad_EFT)*App) &
                      & + a*Av*(-4._dl*(3._dl*a*Omegam_L + 5._dl*Omegarad_EFT)*Ap+ a*((-(a**4*Omegavac_EFT) + 5._dl*a*Omegam_L + 7._dl*Omegarad_EFT)*App &
                      & - a*(a**4*Omegavac_EFT + a*Omegam_L + Omegarad_EFT)*Appp)) +a**2*(2._dl*(3._dl*a*Omegam_L + 5._dl*Omegarad_EFT)*Ap**2 &
                      & + a**2*(a**4*Omegavac_EFT + a*Omegam_L + Omegarad_EFT)*App**2 +a*Ap*((a**4*Omegavac_EFT - 5._dl*a*Omegam_L - 7._dl*Omegarad_EFT)*App &
                      & + a*(a**4*Omegavac_EFT + a*Omegam_L + Omegarad_EFT)*Appp))))/a**6

                vfunc= -Av**2 *(1._dl-eps2)*Omegavac_EFT+Av**2 *(1._dl-eps2)*Omegarad_EFT/(3._dl*a**4) - Av**4 * Omegarad_EFT/(3.*a**4)+Av**2 *(1-eps2) &
                    & * (eps2+2._dl*a*eps2p/(1._dl-eps2))*(Omegam_L/a**3 + Omegarad_EFT/a**4 + Omegavac_EFT)/3._dl

                vfuncp = (Av*(4._dl*Omegarad_EFT*Av**3 - 4._dl*a*Omegarad_EFT*Av**2*Ap + 2._dl*a*Ap*((1._dl - eps2) &
                       &  *(-3._dl*a**4*Omegavac_EFT + Omegarad_EFT + (a**4*Omegavac_EFT + a*Omegam_L + Omegarad_EFT)*eps2) &
                       & + 2._dl*a*(a**4*Omegavac_EFT + a*Omegam_L + Omegarad_EFT)*eps2p) + Av*(-4._dl*Omegarad_EFT &
                       & + (3._dl*a*Omegam_L + 4._dl*Omegarad_EFT)*eps2**2 + a*eps2*(-3._dl*Omegam_L - 2._dl*(a**4*Omegavac_EFT + a*Omegam_L + Omegarad_EFT)*eps2p) &
                       & + a*((6._dl*a**4*Omegavac_EFT - 3._dl*a*Omegam_L - 6._dl*Omegarad_EFT)*eps2p + 2._dl*a*(a**4*Omegavac_EFT + a*Omegam_L + Omegarad_EFT)*eps2pp))))/(3._dl*a**5)

                vfuncpp= -(20._dl*Omegarad_EFT*Av**4 + 4._dl*Omegarad_EFT*Av**2*(-5._dl + 3._dl*a**2 * Ap**2) + 4._dl*a*Omegarad_EFT*Av**3 * (-8._dl*Ap + a*App) &
                    & + a**2* (2._dl*(9._dl*a*Omegam_L + 14._dl*Omegarad_EFT)*Ap**2 + 4._dl*a**2*(a**4 *Omegavac_EFT + a*Omegam_L + Omegarad_EFT)*App**2 &
                    & + a*Ap*((-2._dl* a**4 *Omegavac_EFT - 17._dl*a*Omegam_L - 22._dl*Omegarad_EFT)*App + 2._dl*a*(a**4 *Omegavac_EFT + a*Omegam_L + Omegarad_EFT)*Appp)) &
                    & + a*Av*(-2._dl*(9._dl*a*Omegam_L + 4._dl*Omegarad_EFT)*Ap + a*(-2._dl*(5._dl*a**4 *Omegavac_EFT - 4._dl*a*Omegam_L + Omegarad_EFT)*App &
                    & + a*((-14._dl*a**4 *Omegavac_EFT + a*Omegam_L + 6._dl*Omegarad_EFT)*Appp - 2._dl*a*(a**4 *Omegavac_EFT + a*Omegam_L + Omegarad_EFT)*Apppp))))/(3._dl*a**6)
                
                ! 2) Get the derivatives of chi w.r.t. x:
                ydot(1)= 2._dl*vfuncp/(ufunc+vfunc) * y(1)*a

                chi= y(1)

                chip= ydot(1)/a
              
                chipp= (2._dl*chi*(-ufuncp *vfuncp+vfuncp**2 +(ufunc + vfunc)*vfuncpp))/(ufunc + vfunc)**2

                K= vfunc/OmegaPhi_0

                Kp= vfuncp/OmegaPhi_0

                dKdchi= (ufunc+vfunc)/(2._dl*chi*OmegaPhi_0)

                d2Kdchi2=  (ufunc+vfunc)*(ufuncp-vfuncp)/(4._dl*chi**2 *OmegaPhi_0*vfuncp)

                d3Kdchi3= ((ufunc + vfunc)*(vfuncp*(ufuncp**2-4._dl*ufuncp*vfuncp+ 3._dl*vfuncp**2 + (ufunc + vfunc)*ufuncpp)-(ufunc + vfunc)*ufuncp* vfuncpp)) &
                          & /(8._dl*OmegaPhi_0*chi**3*vfuncp**3)

                adotoa= a *params_cache%h0_Mpc *(Av/(1._dl-eps2))*sqrt(Omegam_EFT/a**3 + Omegarad_EFT/a**4 + OmegaPhi_0* Av**(-4)*(2._dl*chi*dKdchi-K))

                ydot(2) = - sqrt( 2._dl*chi* Mnorm*grhom)* a/(Av*adotoa) 

                Hdot= -((params_cache%h0_Mpc)**2*(-2._dl*a*(a*Omegam_EFT + Omegarad_EFT)*Av**4 *(-1._dl + eps2)*Ap + 2._dl*a**5* OmegaPhi_0*(-1._dl + eps2)*(2._dl*chi*dKdchi - K) &                   
                    & * Ap + Av**5 * ((a*Omegam_EFT + 2._dl *Omegarad_EFT)*(-1._dl + eps2) + 2._dl*a*(a*Omegam_EFT + Omegarad_EFT)*eps2p) & 
                    & + a**4 * OmegaPhi_0 *Av*(-2._dl*(-1._dl + eps2)*(-K + a*dKdchi*chip + chi*(2._dl*dKdchi + a*chip*d2Kdchi2)) + 2._dl*a*(2._dl*chi*dKdchi - K)*eps2p & 
                    & + a*(-1._dl+ eps2)*Kp)))/(2._dl*a**2 *Av**3 *(-1._dl + eps2)**3)

                adotoaL = a * params_cache%h0_Mpc * sqrt(Omegam_EFT/a**3 + Omegarad_EFT/a**4 + 1._dl-Omegam_EFT - Omegarad_EFT)
     
                HdotL = -1.5_dl*a**2 * params_cache%h0_Mpc**2 * (Omegam_EFT/a**3 + 4._dl*Omegarad_EFT/(3._dl*a**4)) + adotoaL**2


            end if 

            rho_phi= Mnorm*grhom*Av**(-4)*(2._dl*chi*dKdchi-K)

            p_phi= Mnorm*grhom*Av**(-4)*K

            eps1= eps2**2 *Av**2 *adotoa**2/( Mnorm*grhom* chi*dKdchi)

            OmegaPhi= rho_phi* OmegaPhi_0*(params_cache%h0_Mpc*Av*a/adotoa)**2 /(Mnorm*grhom)

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
        subroutine output( num_eq, ind, x, y )

            implicit none

            integer , intent(in)                     :: num_eq !< number of equations in the ODE system.
            integer , intent(in)                     :: ind    !< index of the EFT functions interpolation tables to fill.
            real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed.
            real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system.

            logical :: is_open

            ! 1) call derivs to make sure everything is initialized at the correct time step:
            call derivs( num_eq, x, y, ydot )

            ! 2) compute the EFT functions:
            ! NB: The Planck mass is included in grhom
            self%EFTOmega%y(ind)     = Av**(-2) - 1._dl
            self%EFTOmega%yp(ind)    = - 2._dl * Av**(-3) *Ap
            self%EFTOmega%ypp(ind)   = 6._dl * Av**(-4) *Ap**2 - 2._dl* Av**(-3)*App
            self%EFTOmega%yppp(ind)  = - 24._dl * Av**(-5) *Ap**3 + 18._dl* Av**(-4)*Ap*App - 2._dl* Av**(-3)*Appp  

            self%EFTLambda%y(ind)  = (a**2)*grhom*Mnorm*K/(Av**4) - 3._dl*(adotoa**2)*(eps2**2)/(Av**2)
            self%EFTLambda%yp(ind) = adotoa*(- 4._dl*a**3*grhom*Mnorm*Ap*K+a**3*grhom*Mnorm*Av*chip* dKdchi +6._dl*Ap*a* Av**2 &
                & * eps2**2 *adotoa**2+ 6._dl* Av**3 * eps2*adotoa*(eps2*(adotoa-Hdot/adotoa)- adotoa &
                & *(eps2-eps2**2 +(a**2 *App/Av))))/Av**5 

            self%EFTc%y(ind)  = (a**2)*grhom*Mnorm*chi*dKdchi/(Av**4) - 3._dl*(adotoa**2)*(eps2**2)/(Av**2)
            self%EFTc%yp(ind) =  adotoa*(- 4._dl*a**3*grhom*Mnorm*Ap*chi*dKdchi+a**3 *grhom*Mnorm*Av*chip*(chi*d2Kdchi2+dKdchi) &
                & + 6._dl*a*Ap* Av**2 *eps2**2 *adotoa**2+ 6._dl* Av**3 * eps2*adotoa*(eps2*(adotoa-Hdot/adotoa) &
                & -adotoa*(eps2-eps2**2 +(a**2 *App/Av))))/Av**5

            self%EFTgamma1%y(ind)  = (Av**(-4))*grhom*Mnorm*chi**2 *d2Kdchi2/((params_cache%h0_Mpc)**2)
            self%EFTgamma1%yp(ind) = (-4._dl *Ap/Av + chip * (2._dl/chi + d3Kdchi3/d2Kdchi2)) *(Av**(-4)) &
                & * grhom*Mnorm*chi**2 *d2Kdchi2/((params_cache%h0_Mpc)**2)

            ! write debug quantities:
            if ( DebugEFTCAMB ) then
                inquire( unit=33, opened=is_open )
                if ( is_open ) then
                    write (33,'(20E15.5)') a, adotoa, (adotoa- adotoaL)/adotoaL, (1+eps1)*Av**2 -1.,  OmegaPhi, p_phi/rho_phi, eps1, Av, K, y(2), chi
                end if
                inquire( unit=44, opened=is_open )
                if ( is_open ) then
                    write (44,'(20E15.5)') a, self%EFTOmega%y(ind),self%EFTOmega%yp(ind),self%EFTOmega%ypp(ind),self%EFTOmega%yppp(ind), self%EFTLambda%y(ind) , self%EFTLambda%yp(ind),self%EFTc%y(ind) ,self%EFTc%yp(ind)
                end if
            end if

        end subroutine

    end subroutine EFTCAMBKmouflageSolveBackgroundEquations

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBKmouflageAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Kmouflage)                     :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)  :: m_lim  

        logical :: EFTCAMBKmouflageAdditionalModelStability           !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBKmouflageAdditionalModelStability = .True.

        m_lim=0._dl
        
        if ( self%alphaU < 0.0_dl ) then
            EFTCAMBKmouflageAdditionalModelStability = .False.
            return
        end if

        if ( .not. self%Kmimic .and. self%gammaU < 1.0_dl) then
            EFTCAMBKmouflageAdditionalModelStability = .False.
            return
        end if

        if (self%m <= 1.0_dl ) then
            print*, 'Chose an appropriate value for parameter m'
            EFTCAMBKmouflageAdditionalModelStability = .False.
            return
        end if

        if ( .not. self%Kmimic .and. self%eps2_0 > 0.0_dl ) then
            EFTCAMBKmouflageAdditionalModelStability = .False.
            return
        else if ( self%Kmimic .and.  self%eps2_0 < 0.0_dl) then 
            EFTCAMBKmouflageAdditionalModelStability = .False.
            return
        end if

        if ( .not. self%Kmimic .and. self%gammaA < 0.0_dl ) then
            EFTCAMBKmouflageAdditionalModelStability = .False.
            return
        else if ( self%Kmimic .and.  self%gammaA < self%eps2_0/(1._dl-self%eps2_0) )then 
            EFTCAMBKmouflageAdditionalModelStability = .False.
            return
        end if

    end function EFTCAMBKmouflageAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_Kmouflage_Mod


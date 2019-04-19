! -------------------------------------------------------------------------------------------------
!
!   EFTCAMB
!
!   Developed and implemented by:
!       Bin Hu (hu@lorentz.leidenuniv.nl), Marco Raveri (mraveri@sissa.it)
!       Noemi Frusciante (fruscian@iap.fr), Alessandra Silvestri (silvestri@lorentz.leidenuniv.nl)
!
!
!   Simple program that will sample parameter space to output the result of the stability
!   of the dark energy / Modified Gravity model considered.
!
!   This program requires compiling with camb but actually does not use the parameter file.
!
!   For more informations about the methods contained in this file see
!   the documentation: arXiv:1405.3590
!
! -------------------------------------------------------------------------------------------------

module StabilityFunctions
    use precision
    use CAMB
    implicit none

    type(CAMBparams), target :: P
    real(dl), pointer :: EFTparam
    real(dl), pointer :: EFTparam2
    real(dl), pointer :: EFTparam3

    logical :: testCls = .false.
    logical :: showS8  = .false.

contains

    ! ----------------------------

    function Test1D(param)
        ! Step function for 1D test purposes.
        implicit none

        real(dl) :: Test1D
        real(dl), intent(in) :: param

        if (param>0._dl) then
            Test1D = 1._dl
        else
            Test1D = 0._dl
        end if

        return
    end function Test1D

    ! ----------------------------

    function Test2D(param1, param2)
        ! Step function for 2D test purposes.
        implicit none

        real(dl) :: Test2D
        real(dl), intent(in) :: param1, param2

        if (param1 > 0._dl.and.param2 > 0._dl) then
            Test2D = 1._dl
        else
            Test2D = 0._dl
        end if

        return
    end function Test2D

    ! ----------------------------

    function Test3D(param1, param2, param3)
        ! Step function for 3D test purposes.
        implicit none

        real(dl) :: Test3D
        real(dl), intent(in) :: param1, param2, param3

        if (param1 > 0._dl.and.param2 > 0._dl.and.param3 > 0._dl) then
            Test3D = 1._dl
        else
            Test3D = 0._dl
        end if

        return
    end function Test3D

    ! ----------------------------

    function Stability1D(param)
        ! Binary function that gives the stability varying 1 parameter.
        use EFTinitialization
        implicit none

        real(dl) :: Stability1D
        real(dl), intent(in) :: param
        logical  :: success

        ! 1) Set the cosmological parameter that is varied:
        EFTparam = param
        ! 2) Call the CAMB subroutine to set cosmological parameters:
        call CAMBParams_Set(P)
        success = .true.
        ! 3) Call the EFT initialization module to check the stability:
        call EFTCAMB_initialization(success)
        ! 4) Write the result:
        if (success) then
            Stability1D = 1._dl
            if (testCls) then
                call CAMB_GetResults(P)
                if (showS8) call Transfer_output_sig8(MT)
            end if
        else
            Stability1D = 0._dl
        end if

        return
    end function Stability1D

    ! ----------------------------

    function Stability2D(param1, param2)
        ! Binary function that gives the stability varying 2 parameters.
        use EFTinitialization
        implicit none

        real(dl) :: Stability2D
        real(dl), intent(in) :: param1, param2
        logical  :: success

        ! 1) Set the cosmological parameter that is varied:
        EFTparam  = param1
        EFTparam2 = param2
        ! 2) Call the CAMB subroutine to set cosmological parameters:
        call CAMBParams_Set(P)
        success = .true.
        ! 3) Call the EFT initialization module to check the stability:
        call EFTCAMB_initialization(success)
        ! 4) Write the result:
        if (success) then
            Stability2D = 1._dl
            if (testCls) then
                call CAMB_GetResults(P)
                if (showS8) call Transfer_output_sig8(MT)
            end if
        else
            Stability2D = 0._dl
        end if

        return
    end function Stability2D

     ! ----------------------------

    function Stability3D(param1, param2, param3)
        ! Binary function that gives the stability varying 2 parameters.
        use EFTinitialization
        implicit none

        real(dl) :: Stability3D
        real(dl), intent(in) :: param1, param2, param3
        logical  :: success

        ! 1) Set the cosmological parameter that is varied:
        EFTparam  = param1
        EFTparam2 = param2
        EFTparam3 = param3
        ! 2) Call the CAMB subroutine to set cosmological parameters:
        call CAMBParams_Set(P)
        success = .true.
        ! 3) Call the EFT initialization module to check the stability:
        call EFTCAMB_initialization(success)
        ! 4) Write the result:
        if (success) then
            Stability3D = 1._dl
            if (testCls) then
                call CAMB_GetResults(P)
                if (showS8) call Transfer_output_sig8(MT)
            end if
        else
            Stability3D = 0._dl
        end if

        return
    end function Stability3D

end module StabilityFunctions

! -------------------------------------------------------------------------------------------------

module SamplersSubroutines
    use precision
    implicit none

contains

    ! ----------------------------

    subroutine Grid_Sampling1D(function, IntMin, IntMax, Points, Results)
        ! Subroutine that performs sampling of a 1D function on a grid of points.
        !
        ! function = function to be sampled                                     (in)
        ! IntMin   = minimum of the interval in which the function is sampled   (in)
        ! IntMax   = maximum of the interval in which the function is sampled   (in)
        ! Points   = number of points to sample                                 (in)
        ! Results  = input: an allocatable array, output: allocated array of    (in out)
        !            dimension (2,Points) containing the points of the form
        !            (x,function(x))
        !
        implicit none

        interface
            function function(x)
                use precision
                real(dl) ::  function
                real(dl), intent(in) :: x
            end function function
        end interface

        real(dl), intent(in)  :: IntMin, IntMax
        integer , intent(in)  :: Points
        real(dl), allocatable :: Results(:,:)

        integer  :: i,j
        real(dl) :: x, y

        ! 1) Deallocate and allocate the results array.
        if (allocated(Results)) deallocate(Results)
        allocate(Results(2,Points))
        ! 2) Do the sampling
        do i = 1, Points
            x = IntMin +REAL(i-1)*(IntMax-IntMin)/REAL(Points-1)
            y = function(x)

            Results(1,i) = x
            Results(2,i) = y
        end do

        return
    end subroutine Grid_Sampling1D

    ! ----------------------------

    subroutine Adaptive_Sampling1D(function, IntMin, IntMax, Points, Results, Recursions)
        ! Subroutine that performs an adaptively refined sampling of a 1D function.
        ! Much more efficient than grid sampling.
        !
        ! function   = function to be sampled                                     (in)
        ! IntMin     = minimum of the interval in which the function is sampled   (in)
        ! IntMax     = maximum of the interval in which the function is sampled   (in)
        ! Points     = number of points of the initial grid sampling              (in)
        ! Results    = input: an allocatable array, output: allocated array       (in out)
        !              containing the points of the form (x,function(x))
        ! Recursions = (OPTIONAL) maximum number of points in which               (in)
        !              every sub-interval is adaptively sampled.
        !
        implicit none

        interface
            function function(x)
                use precision
                real(dl) ::  function
                real(dl), intent(in) :: x
            end function function
        end interface

        real(dl), intent(in)  :: IntMin, IntMax
        integer , intent(in)  :: Points
        integer , optional    :: Recursions

        real(dl), allocatable :: Results(:,:)
        real(dl), allocatable :: Temp(:,:), Temp2(:,:)

        integer  :: i,j
        integer  :: GlobalCounter, LocalCounter
        real(dl) :: x, y, x1, x2, y1, y2

        ! 1) If the maximum number of recursions is not in input set it to 5
        if (.not.present(Recursions)) Recursions = 5
        ! 2) Allocate Temp as the maximum number of points possible given Recursions
        if (allocated(Temp))  deallocate(Temp)
        if (allocated(Temp2)) deallocate(Temp2)
        allocate(Temp(2,Points*(Recursions+1)),Temp2(2,Points))
        ! 3) Do the sampling
        GlobalCounter = 1
        do i = 1, Points
            x = IntMin +REAL(i-1)*(IntMax-IntMin)/REAL(Points-1)
            y = function(x)
            Temp2(1,i) = x
            Temp2(2,i) = y
            if (i>1) then
                LocalCounter  = 0
                x1 = Temp2(1,i-1)
                x2 = Temp2(1,i)
                y1 = Temp2(2,i-1)
                y2 = Temp2(2,i)
                call Interval_Recursive_Sampler(function, x1, x2, y1, y2, LocalCounter, Recursions,&
                    & Temp, Points*(Recursions+1), GlobalCounter)
            end if
            Temp(1,GlobalCounter) = x
            Temp(2,GlobalCounter) = y
            GlobalCounter = GlobalCounter + 1
        end do

        ! 4) Write out results
        if (allocated(Results)) deallocate(Results)
        allocate(Results(2,GlobalCounter-1))
        do i = 1, GlobalCounter-1
            Results(1,i) = Temp(1,i)
            Results(2,i) = Temp(2,i)
        end do

        ! 5) Free memory
        if (allocated(Temp))  deallocate(Temp)
        if (allocated(Temp2)) deallocate(Temp2)

        return
    end subroutine Adaptive_Sampling1D

    ! ----------------------------

    recursive subroutine Interval_Recursive_Sampler(function, x1, x2, y1, y2, counter, MaxRecursion,&
        & ResArray, ResArrayLength, ResArrayInd)
        ! Subdivides an intervall and decide wether to sample it recursively.
        ! This works good for our purposes but not so well for a generic function.
        ! The criterium for refining is based on a combination of first and second derivatives.
        implicit none

        external function
        real(dl) :: function

        real(dl), intent(inout) :: x1, x2, y1, y2
        integer , intent(inout) :: counter
        integer , intent(in)    :: MaxRecursion
        integer  :: ResArrayLength, ResArrayInd
        real(dl) :: ResArray(2,ResArrayLength)


        real(dl), parameter :: tol = 0.1_dl
        real(dl), parameter :: cut = 1000._dl

        real(dl) :: xm, ym, ep
        real(dl) :: der2

        if (counter >= MaxRecursion) return

        ep = ABS(y2-y1)/ABS(x2-x1)

        if (ep>tol.and.ep<cut) then
            xm = 0.5_dl*(x1+x2)
            ym = function(xm)
            der2 = (y1 -2._dl*ym +y2)/(x2-xm)**2

            counter = counter + 1

            ResArray(1,ResArrayInd) = xm
            ResArray(2,ResArrayInd) = ym

            ResArrayInd = ResArrayInd +1

            if (ABS(der2)>0.1_dl*tol) then
                call Interval_Recursive_Sampler(function, x1, xm, y1, ym, counter, MaxRecursion,&
                    &ResArray, ResArrayLength, ResArrayInd)
                call Interval_Recursive_Sampler(function, xm, x2, ym, y2, counter, MaxRecursion,&
                    &ResArray, ResArrayLength, ResArrayInd)
            end if
        else
            return
        end if

    end subroutine Interval_Recursive_Sampler

    ! ----------------------------

    subroutine Grid_Sampling2D(function, IntMin1, IntMax1, IntMin2, IntMax2, Points, Results)
        ! function = function to be sampled                                     (in)
        ! IntMin1  = minimum of the interval in which the function is sampled   (in)
        !            in the first parameter.
        ! IntMax1  = maximum of the interval in which the function is sampled   (in)
        !            in the first parameter.
        ! IntMin2  = minimum of the interval in which the function is sampled   (in)
        !            in the second parameter.
        ! IntMax2  = maximum of the interval in which the function is sampled   (in)
        !            in the second parameter.
        ! Points   = number of points in which both intervals are subdivided    (in)
        !            the total number of sampled points would be Points**2
        ! Results  = input: an allocatable array, output: allocated array of    (in out)
        !            dimension (2,Points) containing the points of the form
        !            (x,function(x))
        implicit none
        ! MRTODO: write the interface for f
        external function
        real(dl) :: function

        real(dl), intent(in)  :: IntMin1, IntMax1, IntMin2, IntMax2
        integer , intent(in)  :: Points
        real(dl), allocatable :: Results(:,:)

        integer  :: i,j,k
        real(dl) :: x, y, z

        ! 1) Deallocate and allocate the results array.
        if (allocated(Results)) deallocate(Results)
        allocate(Results(3,Points**2))
        ! 2) Do the sampling
        k=0
        do i = 1, Points
            do j = 1, Points

                x = IntMin1 +REAL(i-1)*(IntMax1-IntMin1)/REAL(Points-1)
                y = IntMin2 +REAL(j-1)*(IntMax2-IntMin2)/REAL(Points-1)
                z = function(x,y)

                k = k +1

                Results(1,k) = x
                Results(2,k) = y
                Results(3,k) = z

            end do
        end do

        return
    end subroutine Grid_Sampling2D

    ! ----------------------------

    subroutine Grid_Sampling3D(function, IntMin1, IntMax1, IntMin2, IntMax2, IntMin3, IntMax3, Points, Results)
        ! function = function to be sampled                                     (in)
        ! IntMin1  = minimum of the interval in which the function is sampled   (in)
        !            in the first parameter.
        ! IntMax1  = maximum of the interval in which the function is sampled   (in)
        !            in the first parameter.
        ! IntMin2  = minimum of the interval in which the function is sampled   (in)
        !            in the second parameter.
        ! IntMax2  = maximum of the interval in which the function is sampled   (in)
        !            in the second parameter.
        ! IntMin3  = minimum of the interval in which the function is sampled   (in)
        !            in the third parameter.
        ! IntMax3  = maximum of the interval in which the function is sampled   (in)
        !            in the third parameter.
        ! Points   = number of points in which both intervals are subdivided    (in)
        !            the total number of sampled points would be Points**2
        ! Results  = input: an allocatable array, output: allocated array of    (in out)
        !            dimension (2,Points) containing the points of the form
        !            (x,function(x))
        implicit none
        ! MRTODO: write the interface for f
        external function
        real(dl) :: function

        real(dl), intent(in)  :: IntMin1, IntMax1, IntMin2, IntMax2, IntMin3, IntMax3
        integer , intent(in)  :: Points
        real(dl), allocatable :: Results(:,:)

        integer  :: i,j,k,m
        real(dl) :: x, y, z, f

        ! 1) Deallocate and allocate the results array.
        if (allocated(Results)) deallocate(Results)
        allocate(Results(4,Points**3))
        ! 2) Do the sampling
        k=0
        do i = 1, Points
            do j = 1, Points
                do m = 1, Points

                    x = IntMin1 +REAL(i-1)*(IntMax1-IntMin1)/REAL(Points-1)
                    y = IntMin2 +REAL(j-1)*(IntMax2-IntMin2)/REAL(Points-1)
                    z = IntMin3 +REAL(m-1)*(IntMax3-IntMin3)/REAL(Points-1)
                    f = function(x,y,z)

                    k = k +1

                    Results(1,k) = x
                    Results(2,k) = y
                    Results(3,k) = z
                    Results(4,k) = f

                end do
            end do
        end do

        return
    end subroutine Grid_Sampling3D

end module SamplersSubroutines

! -------------------------------------------------------------------------------------------------

module Systematic_Model_Sampling
    use CAMB
    use EFTinitialization
    use StabilityFunctions
    use SamplersSubroutines

    use omp_lib

    implicit none

contains

    ! ----------------------------

    subroutine EFT_Clean_Params(P)
        ! Subroutines that resets the values of the EFT parameters
        implicit none
        type(CAMBparams) :: P
        ! 1) Clear flags:
        P%EFTflag = 0
        P%EFTwDE  = 0
        P%PureEFTmodelOmega  = 0
        P%PureEFTmodelGamma1 = 0
        P%PureEFTmodelGamma2 = 0
        P%PureEFTmodelGamma3 = 0
        P%PureEFTmodelGamma4 = 0
        P%PureEFTmodelGamma5 = 0
        P%PureEFTmodelGamma6 = 0
        P%DesignerEFTmodel   = 0
        P%FullMappingEFTmodel   = 0
        P%HoravaSolarSystem  = .false.
        P%PureEFTHorndeski   = .false.
        P%RPHmassPmodel      = 0
        P%RPHkineticitymodel = 0
        P%RPHbraidingmodel   = 0
        P%RPHtensormodel     = 0
        ! 2) Clear model parameters:
        P%EFTw0 = 0._dl
        P%EFTwa = 0._dl
        P%EFTwn = 0._dl
        P%EFTwat= 0._dl
        P%EFtw2 = 0._dl
        P%EFTw3 = 0._dl
        P%EFTOmega0    = 0._dl
        P%EFTOmegaExp  = 0._dl
        P%EFTGamma10   = 0._dl
        P%EFTGamma1Exp = 0._dl
        P%EFTGamma20   = 0._dl
        P%EFTGamma2Exp = 0._dl
        P%EFTGamma30   = 0._dl
        P%EFTGamma3Exp = 0._dl
        P%EFTGamma40   = 0._dl
        P%EFTGamma4Exp = 0._dl
        P%EFTGamma50   = 0._dl
        P%EFTGamma5Exp = 0._dl
        P%EFTGamma60   = 0._dl
        P%EFTGamma6Exp = 0._dl
        P%EFTB0        = 0._dl
        P%RPHmassP0        = 0._dl
        P%RPHmassPexp      = 0._dl
        P%RPHkineticity0   = 0._dl
        P%RPHkineticityexp = 0._dl
        P%RPHbraiding0     = 0._dl
        P%RPHbraidingexp   = 0._dl
        P%RPHtensor0       = 0._dl
        P%RPHtensorexp     = 0._dl
        P%Horava_xi        = 0._dl
        P%Horava_lambda    = 0._dl
        P%Horava_eta       = 0._dl

    end subroutine EFT_Clean_Params

    ! ----------------------------

    subroutine SamplePureEFT1D(NumPoints, IntMin, IntMax, OutRoot, MaxRecursions)
        ! Subroutine to test at once every 1 parameter, pure EFT, extension of the standard model.
        implicit none

        integer , intent(in) :: NumPoints
        real(dl), intent(in) :: IntMin, IntMax

        integer, optional :: MaxRecursions

        character(LEN=100), intent(in) :: OutRoot
        character(LEN=100)             :: OutDir
        real(dl), allocatable :: res(:,:)

        ! timings:
        real(dl) :: t1, t2, total_time = 0._dl
        integer  :: total_number = 0

        ! 1) Set MaxRecursions
        if (.not.present(MaxRecursions)) MaxRecursions = 5
        ! 2) Set the output directory
        OutDir = TRIM(OutRoot)//'1_PureEFT_1P/'
        call system('mkdir '//OutDir)
        ! 3) Some feedback
        write(*,*) 'Starting Pure EFT one parameter stability check:'
        write(*,*) ' results files placed in: '//TRIM(OutDir)

        ! 1) Pure EFT Omega, constant
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelOmega  = 1
        EFTparam => P%EFTOmega0
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Omega_constant.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Omega constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time   = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 2) Pure EFT Omega, linear
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelOmega  = 2
        EFTparam => P%EFTOmega0
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Omega_linear.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Omega linear in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 3) Pure EFT Gamma1, constant
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma1 = 1
        EFTparam => P%EFTGamma10
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma1_constant.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Gamma1 constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

         ! 4) Pure EFT Gamma1, linear
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma1 = 2
        EFTparam => P%EFTGamma10
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma1_linear.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Gamma1 linear in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 5) Pure EFT Gamma2, constant
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma2 = 1
        EFTparam => P%EFTGamma20
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma2_constant.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Gamma2 constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

         ! 6) Pure EFT Gamma2, linear
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma2 = 2
        EFTparam => P%EFTGamma20
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma2_linear.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Gamma2 linear in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 7) Pure EFT Gamma3, constant
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma3 = 1
        EFTparam => P%EFTGamma30
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma3_constant.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Gamma3 constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

         ! 8) Pure EFT Gamma3, linear
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma3 = 2
        EFTparam => P%EFTGamma30
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma3_linear.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Gamma3 linear in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 9) Pure EFT Gamma4, constant
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma4 = 1
        EFTparam => P%EFTGamma40
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma4_constant.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Gamma4 constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

         ! 10) Pure EFT Gamma4, linear
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma4 = 2
        EFTparam => P%EFTGamma40
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma4_linear.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Gamma4 linear in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 11) Pure EFT Gamma5, constant
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma5 = 1
        EFTparam => P%EFTGamma50
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma5_constant.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Gamma5 constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

         ! 12) Pure EFT Gamma5, linear
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma5 = 2
        EFTparam => P%EFTGamma50
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma5_linear.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Gamma5 linear in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 13) Pure EFT Gamma6, constant
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma6 = 1
        EFTparam => P%EFTGamma60
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma6_constant.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Gamma6 constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

         ! 14) Pure EFT Gamma6, linear
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma6 = 2
        EFTparam => P%EFTGamma60
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma6_linear.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') '  Done EFT Gamma6 linear in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        write(*,*) 'Pure EFT 1D stability check done.'
        write(*,'(a,F9.2,a,F9.5,a)') ' Total time: ', total_time, &
                                       &' sec. Average time per call: ', total_time/total_number, ' sec.'

    end subroutine SamplePureEFT1D

    ! ----------------------------

    subroutine SamplePureEFTwCDM_2D(NumPoints, IntMin1, IntMax1, IntMin2, IntMax2, OutRoot)
        ! Subroutine to test at once every 1 parameter, pure EFT, extension of the standard model
        ! on wCDM background.
        implicit none

        integer , intent(in) :: NumPoints
        real(dl), intent(in) :: IntMin1, IntMax1, IntMin2, IntMax2
        character(LEN=100), intent(in) :: OutRoot
        character(LEN=100)             :: OutDir

        real(dl), allocatable :: res(:,:)

        ! timings:
        real(dl) :: t1, t2, total_time = 0._dl
        integer  :: total_number = 0

        ! 2) Set the output directory
        OutDir = TRIM(OutRoot)//'2_PureEFT_wCDM_2D/'
        call system('mkdir '//OutDir)
        ! 3) Some feedback
        write(*,*) 'Starting Pure EFT one parameter on wCDM background stability check:'
        write(*,*) ' results files placed in: '//TRIM(OutDir)

        ! 1) Pure EFT Omega, constant, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelOmega  = 1
        EFTparam  => P%EFTOmega0
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Omega_const_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Omega const on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 2) Pure EFT Omega, linear, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelOmega  = 2
        EFTparam  => P%EFTOmega0
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Omega_linear_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Omega linear on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 3) Pure EFT Gamma1, constant, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelGamma1 = 1
        EFTparam  => P%EFTGamma10
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma1_const_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma1 constant on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 4) Pure EFT Gamma1, linear, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelGamma1 = 2
        EFTparam  => P%EFTGamma10
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma1_linear_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma1 linear on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 5) Pure EFT Gamma2, constant, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelGamma2 = 1
        EFTparam  => P%EFTGamma20
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma2_const_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma2 const on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 6) Pure EFT Gamma2, linear, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelGamma2 = 2
        EFTparam  => P%EFTGamma20
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma2_linear_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma2 linear on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

         ! 7) Pure EFT Gamma3, constant, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelGamma3 = 1
        EFTparam  => P%EFTGamma30
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma3_const_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma3 const on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 8) Pure EFT Gamma3, linear, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelGamma3 = 2
        EFTparam  => P%EFTGamma30
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma3_linear_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma3 linear on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 9) Pure EFT Gamma4, constant, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelGamma4 = 1
        EFTparam  => P%EFTGamma40
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma4_const_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma4 const on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 10) Pure EFT Gamma4, linear, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelGamma4 = 2
        EFTparam  => P%EFTGamma40
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma4_linear_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma4 linear on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 11) Pure EFT Gamma5, constant, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelGamma5 = 1
        EFTparam  => P%EFTGamma50
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma5_const_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma5 const on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 12) Pure EFT Gamma5, linear, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelGamma5 = 2
        EFTparam  => P%EFTGamma50
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma5_linear_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma5 linear on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 13) Pure EFT Gamma6, constant, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelGamma6 = 1
        EFTparam  => P%EFTGamma60
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma6_const_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma6 const on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 14) Pure EFT Gamma6, linear, wCDM
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 1
        P%PureEFTmodelGamma6 = 2
        EFTparam  => P%EFTGamma60
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma6_linear_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma6 linear on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        write(*,*) 'Pure EFT stability on wCDM check done.'
        write(*,'(a,F9.2,a,F9.5,a)') ' Total time: ', total_time, &
                                       &' sec. Average time per call: ', total_time/total_number, ' sec.'

    end subroutine SamplePureEFTwCDM_2D

    ! ----------------------------

    subroutine SamplePureEFT_PowerLaw(NumPoints, IntMin1, IntMax1, IntMin2, IntMax2, OutRoot)
        ! Subroutine to test at once every pure EFT power law model.
        implicit none

        integer , intent(in) :: NumPoints
        real(dl), intent(in) :: IntMin1, IntMax1, IntMin2, IntMax2
        character(LEN=100), intent(in) :: OutRoot
        character(LEN=100)             :: OutDir

        real(dl), allocatable :: res(:,:)

        ! timings:
        real(dl) :: t1, t2, total_time = 0._dl
        integer  :: total_number = 0

        ! 2) Set the output directory
        OutDir = TRIM(OutRoot)//'3_PureEFT_PowerLaw/'
        call system('mkdir '//OutDir)
        ! 3) Some feedback
        write(*,*) 'Starting Pure EFT power law stability check:'
        write(*,*) ' results files placed in: '//TRIM(OutDir)

        ! 1) Pure EFT Omega power law
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelOmega  = 3
        EFTparam  => P%EFTOmega0
        EFTparam2 => P%EFTOmegaExp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Omega_PL.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Omega power law in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 2) Pure EFT Gamma1 power law
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelGamma1 = 3
        EFTparam  => P%EFTGamma10
        EFTparam2 => P%EFTGamma1Exp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma1_PL.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma1 power law in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 3) Pure EFT Gamma2 power law
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelGamma2 = 3
        EFTparam  => P%EFTGamma20
        EFTparam2 => P%EFTGamma2Exp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma2_PL.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma2 power law in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 4) Pure EFT Gamma3 power law
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelGamma3 = 3
        EFTparam  => P%EFTGamma30
        EFTparam2 => P%EFTGamma3Exp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma3_PL.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma3 power law in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 5) Pure EFT Gamma4 power law
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelGamma4 = 3
        EFTparam  => P%EFTGamma40
        EFTparam2 => P%EFTGamma4Exp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma4_PL.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma4 power law in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 6) Pure EFT Gamma5 power law
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelGamma5 = 3
        EFTparam  => P%EFTGamma50
        EFTparam2 => P%EFTGamma5Exp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma5_PL.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma5 power law in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 7) Pure EFT Gamma6 power law
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelGamma6 = 3
        EFTparam  => P%EFTGamma60
        EFTparam2 => P%EFTGamma6Exp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma6_PL.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma6 power law in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        write(*,*) 'Pure EFT power law check done.'
        write(*,'(a,F9.2,a,F9.5,a)') ' Total time: ', total_time, &
                                       &' sec. Average time per call: ', total_time/total_number, ' sec.'

    end subroutine SamplePureEFT_PowerLaw

    ! ----------------------------

    subroutine SamplePureEFT_Exponential(NumPoints, IntMin1, IntMax1, IntMin2, IntMax2, OutRoot)
        ! Subroutine to test at once every pure EFT exponential model.
        implicit none

        integer , intent(in) :: NumPoints
        real(dl), intent(in) :: IntMin1, IntMax1, IntMin2, IntMax2
        character(LEN=100), intent(in) :: OutRoot
        character(LEN=100)             :: OutDir

        real(dl), allocatable :: res(:,:)

        ! timings:
        real(dl) :: t1, t2, total_time = 0._dl
        integer  :: total_number = 0

        ! 1) Set the output directory
        OutDir = TRIM(OutRoot)//'4_PureEFT_Exponential/'
        call system('mkdir '//OutDir)
        ! 2) Some feedback
        write(*,*) 'Starting Pure EFT exponential models stability check:'
        write(*,*) ' results files placed in: '//TRIM(OutDir)

        ! 1) Pure EFT Omega exponential
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelOmega  = 4
        EFTparam  => P%EFTOmega0
        EFTparam2 => P%EFTOmegaExp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Omega_exp.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Omega exponential in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 2) Pure EFT Gamma1 exponential
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelGamma1 = 4
        EFTparam  => P%EFTGamma10
        EFTparam2 => P%EFTGamma1Exp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma1_exp.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma1 exponential in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 3) Pure EFT Gamma2 exponential
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelGamma2 = 4
        EFTparam  => P%EFTGamma20
        EFTparam2 => P%EFTGamma2Exp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma2_exp.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma2 exponential in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 4) Pure EFT Gamma3 exponential
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelGamma3 = 4
        EFTparam  => P%EFTGamma30
        EFTparam2 => P%EFTGamma3Exp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma3_exp.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma3 exponential in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 5) Pure EFT Gamma4 exponential
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelGamma4 = 4
        EFTparam  => P%EFTGamma40
        EFTparam2 => P%EFTGamma4Exp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma4_exp.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma4 exponential in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 6) Pure EFT Gamma5 exponential
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelGamma5 = 4
        EFTparam  => P%EFTGamma50
        EFTparam2 => P%EFTGamma5Exp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma5_exp.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma5 exponential in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 7) Pure EFT Gamma6 exponential
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%EFTwDE  = 0
        P%PureEFTmodelGamma6 = 4
        EFTparam  => P%EFTGamma60
        EFTparam2 => P%EFTGamma6Exp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma6_exp.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Gamma6 exponential in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        write(*,*) 'Pure EFT exponential check done.'
        write(*,'(a,F9.2,a,F9.5,a)') ' Total time: ', total_time, &
                                       &' sec. Average time per call: ', total_time/total_number, ' sec.'

    end subroutine SamplePureEFT_Exponential

    ! ----------------------------

    subroutine SampleDesigner1D(NumPoints, OutRoot, MaxRecursions)
        ! subroutine that test the stability of 1D designer models.
        implicit none

        integer , intent(in) :: NumPoints
        integer, optional    :: MaxRecursions

        character(LEN=100), intent(in) :: OutRoot
        character(LEN=100)             :: OutDir

        real(dl), allocatable :: res(:,:)

        ! timings:
        real(dl) :: t1, t2, total_time = 0._dl
        integer  :: total_number = 0

        ! 1) Set MaxRecursions
        if (.not.present(MaxRecursions)) MaxRecursions = 5
        ! 2) Set the output directory
        OutDir = TRIM(OutRoot)//'5_Designer_1D/'
        call system('mkdir '//OutDir)
        ! 3) Some feedback
        write(*,*) 'Starting Pure EFT one parameter stability check:'
        write(*,*) ' results files placed in: '//TRIM(OutDir)

        ! 1) designer f(R) on LCDM:
        call EFT_Clean_Params(P)
        P%EFTflag = 2
        P%DesignerEFTmodel   = 1
        EFTparam => P%EFTB0
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, -2._dl, 2._dl, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        !call Grid_Sampling1D(Stability1D, -2._dl, 2._dl, NumPoints, res)
        call CreateTxtFile(TRIM(OutDir)//'1_des_fR_LCDM.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done designer f(R) on LCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 2) designer mc quintessence on wCDM:
        call EFT_Clean_Params(P)
        P%EFTflag = 2
        P%EFTwDE  = 1
        P%DesignerEFTmodel   = 2
        EFTparam => P%EFTw0
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, -2._dl, 0._dl, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        !call Grid_Sampling1D(Stability1D, -2._dl, 2._dl, NumPoints, res)
        call CreateTxtFile(TRIM(OutDir)//'1_des_5e_wCDM.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done designer mc quintessence on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        write(*,*) 'Designer 1D check done.'
        write(*,'(a,F9.2,a,F9.5,a)') ' Total time: ', total_time, &
                                       &' sec. Average time per call: ', total_time/total_number, ' sec.'

    end subroutine SampleDesigner1D

    ! ----------------------------

    subroutine SampleDesigner2D(NumPoints, OutRoot)
        ! subroutine that test the stability of 2D designer models.
        implicit none

        integer , intent(in) :: NumPoints

        character(LEN=100), intent(in) :: OutRoot
        character(LEN=100)             :: OutDir

        real(dl), allocatable :: res(:,:)

        ! timings:
        real(dl) :: t1, t2, total_time = 0._dl
        integer  :: total_number = 0

        ! 1) Set the output directory
        OutDir = TRIM(OutRoot)//'6_Designer_2D/'
        call system('mkdir '//OutDir)
        ! 2) Some feedback
        write(*,*) 'Starting 2 params designer models stability check:'
        write(*,*) ' results files placed in: '//TRIM(OutDir)

        ! 1) designer f(R) on wCDM:
        call EFT_Clean_Params(P)
        P%EFTflag = 2
        P%EFTwDE  = 1
        P%DesignerEFTmodel   = 1
        EFTparam  => P%EFTB0
        EFTparam2 => P%EFTw0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, -2._dl, 2._dl, -1.5_dl,0._dl, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_des_fR_wCDM.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done designer f(R) on wCDM in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 2) designer mc quintessence on CPL:
        call EFT_Clean_Params(P)
        P%EFTflag = 2
        P%EFTwDE  = 2
        P%DesignerEFTmodel   = 2
        EFTparam  => P%EFTw0
        EFTparam2 => P%EFTwa
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, -2._dl, 0._dl, -1._dl, 1._dl, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_des_5e_CPL.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done designer mc quintessence on CPL in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        write(*,*) 'Designer 2D check done.'
        write(*,'(a,F9.2,a,F9.5,a)') ' Total time: ', total_time, &
                                       &' sec. Average time per call: ', total_time/total_number, ' sec.'

    end subroutine SampleDesigner2D

    ! ----------------------------

    subroutine SampleRPH1D(NumPoints, IntMin, IntMax, OutRoot, MaxRecursions)
        ! Subroutine to test at once every 1 parameter RPH extension of the standard model.
        implicit none

        integer , intent(in) :: NumPoints
        real(dl), intent(in) :: IntMin, IntMax

        integer, optional :: MaxRecursions

        character(LEN=100), intent(in) :: OutRoot
        character(LEN=100)             :: OutDir
        real(dl), allocatable :: res(:,:)

        ! timings:
        real(dl) :: t1, t2, total_time = 0._dl
        integer  :: total_number = 0

        ! 1) Set MaxRecursions
        if (.not.present(MaxRecursions)) MaxRecursions = 5
        ! 2) Set the output directory
        OutDir = TRIM(OutRoot)//'7_RPH_1P/'
        call system('mkdir '//OutDir)
        ! 3) Some feedback
        write(*,*) 'Starting RPH one parameter stability check:'
        write(*,*) ' results files placed in: '//TRIM(OutDir)

        ! 1) RPH Planck Mass constant:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHmassPmodel  = 1
        EFTparam => P%RPHmassP0
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_RPH_MassP_const.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done RPH MassP constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 2) RPH kineticity constant:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHkineticitymodel = 1
        EFTparam => P%RPHkineticity0
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_RPH_Kineticity_const.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done RPH kineticity constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 3) RPH braiding constant:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHbraidingmodel = 1
        EFTparam => P%RPHbraiding0
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_RPH_Braiding_const.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done RPH braiding constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

         ! 4) RPH tensor constant:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHtensormodel = 1
        EFTparam => P%RPHtensor0
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_RPH_Tensor_const.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done RPH tensor constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        write(*,*) 'RPH 1D stability check done.'
        write(*,'(a,F9.2,a,F9.5,a)') ' Total time: ', total_time, &
                                       &' sec. Average time per call: ', total_time/total_number, ' sec.'

    end subroutine SampleRPH1D

    ! ----------------------------

    subroutine SampleRPH2D(NumPoints, IntMin1, IntMax1, IntMin2, IntMax2, OutRoot)
        ! Subroutine to test at once every 2 parameter RPH extension of the standard model.
        implicit none

        integer , intent(in) :: NumPoints
        real(dl), intent(in) :: IntMin1, IntMax1, IntMin2, IntMax2
        character(LEN=100), intent(in) :: OutRoot
        character(LEN=100)             :: OutDir

        real(dl), allocatable :: res(:,:)

        ! timings:
        real(dl) :: t1, t2, total_time = 0._dl
        integer  :: total_number = 0

        ! 1) Set the output directory
        OutDir = TRIM(OutRoot)//'8_RPH_2P/'
        call system('mkdir '//OutDir)
        ! 2) Some feedback
        write(*,*) 'Starting RPH two parameter:'
        write(*,*) ' results files placed in: '//TRIM(OutDir)

        ! 1) RPH Planck Mass power law:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHmassPmodel    = 2
        EFTparam  => P%RPHmassP0
        EFTparam2 => P%RPHmassPexp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_MassP_PL.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done MassP power law in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 2) RPH kineticity power law:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHkineticitymodel    = 2
        EFTparam  => P%RPHkineticity0
        EFTparam2 => P%RPHkineticityexp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Kineticity_PL.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Kineticity power law in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 3) RPH braiding power law:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHbraidingmodel    = 2
        EFTparam  => P%RPHbraiding0
        EFTparam2 => P%RPHbraidingexp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Braiding_PL.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Braiding power law in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 4) RPH tensor power law:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHtensormodel    = 2
        EFTparam  => P%RPHtensor0
        EFTparam2 => P%RPHtensorexp
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Tensor_PL.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Tensor power law in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        write(*,*) 'RPH power law stability check done.'
        write(*,'(a,F9.2,a,F9.5,a)') ' Total time: ', total_time, &
                                       &' sec. Average time per call: ', total_time/total_number, ' sec.'

    end subroutine SampleRPH2D

    ! ----------------------------

    subroutine SampleRPH2D_constant(NumPoints, IntMin1, IntMax1, IntMin2, IntMax2, OutRoot)
        ! Subroutine to test at once every combination of constant RPH functions.
        implicit none

        integer , intent(in) :: NumPoints
        real(dl), intent(in) :: IntMin1, IntMax1, IntMin2, IntMax2
        character(LEN=100), intent(in) :: OutRoot
        character(LEN=100)             :: OutDir

        real(dl), allocatable :: res(:,:)

        ! timings:
        real(dl) :: t1, t2, total_time = 0._dl
        integer  :: total_number = 0

        ! 1) Set the output directory
        OutDir = TRIM(OutRoot)//'9_RPH_2P_const/'
        call system('mkdir '//OutDir)
        ! 2) Some feedback
        write(*,*) 'Starting RPH two parameter constant:'
        write(*,*) ' results files placed in: '//TRIM(OutDir)

        ! 1) RPH Planck Mass and kineticity:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHmassPmodel      = 1
        P%RPHkineticitymodel = 1
        EFTparam  => P%RPHmassP0
        EFTparam2 => P%RPHkineticity0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_MassP_Kinet.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done MassP and Kineticity constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 2) RPH Planck Mass and braiding:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHmassPmodel      = 1
        P%RPHbraidingmodel   = 1
        EFTparam  => P%RPHmassP0
        EFTparam2 => P%RPHbraiding0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_MassP_Braiding.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done MassP and Braiding constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 3) RPH Planck Mass and tensor:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHmassPmodel      = 1
        P%RPHtensormodel     = 1
        EFTparam  => P%RPHmassP0
        EFTparam2 => P%RPHtensor0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_MassP_Tensor.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done MassP and Tensor constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 4) RPH kineticity and braiding:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHkineticitymodel = 1
        P%RPHbraidingmodel  = 1
        EFTparam  => P%RPHkineticity0
        EFTparam2 => P%RPHbraiding0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Kinet_Braiding.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Kineticity and Braiding constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 5) RPH kineticity and tensor:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHkineticitymodel = 1
        P%RPHtensormodel     = 1
        EFTparam  => P%RPHkineticity0
        EFTparam2 => P%RPHtensor0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Kinet_Tensor.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Kineticity and Tensor constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 6) RPH Tensor and Braiding:
        call EFT_Clean_Params(P)
        P%EFTflag = 3
        P%AltParEFTmodel = 1
        P%RPHtensormodel     = 1
        P%RPHbraidingmodel  = 1
        EFTparam  => P%RPHtensor0
        EFTparam2 => P%RPHbraiding0
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Tensor_Braiding.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Tensor and Braiding constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        write(*,*) 'RPH constant combinations stability check done.'
        write(*,'(a,F9.2,a,F9.5,a)') ' Total time: ', total_time, &
                                       &' sec. Average time per call: ', total_time/total_number, ' sec.'

    end subroutine SampleRPH2D_constant

    ! ----------------------------

    subroutine SamplePureEFT_Hordenski_1D(NumPoints, IntMin, IntMax, OutRoot, MaxRecursions)
        ! Subroutine to test at once every 1 parameter RPH extension of the standard model.
        implicit none

        integer , intent(in) :: NumPoints
        real(dl), intent(in) :: IntMin, IntMax

        integer, optional :: MaxRecursions

        character(LEN=100), intent(in) :: OutRoot
        character(LEN=100)             :: OutDir
        real(dl), allocatable :: res(:,:)

        ! timings:
        real(dl) :: t1, t2, total_time = 0._dl
        integer  :: total_number = 0

        ! 1) Set MaxRecursions
        if (.not.present(MaxRecursions)) MaxRecursions = 5
        ! 2) Set the output directory
        OutDir = TRIM(OutRoot)//'10_PEFT_Horn_1D/'
        call system('mkdir '//OutDir)
        ! 3) Some feedback
        write(*,*) 'Starting Pure EFT Horndeski one parameter stability check:'
        write(*,*) ' results files placed in: '//TRIM(OutDir)

        ! Notice that nothing changes with respect to Pure EFT 1D but the fact that in PEFT Horndeski
        ! there is no choice for Gamma4, Gamma5 and Gamma6 that are related to Gamma3.

        ! 1) Pure EFT Gamma3 Horndeski, constant
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma3 = 1
        P%PureEFTHorndeski   = .true.
        EFTparam => P%EFTGamma30
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma3_Hor_constant.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done EFT Gamma3 Horndeski constant in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

         ! 8) Pure EFT Gamma3 Horndeski, linear
        call EFT_Clean_Params(P)
        P%EFTflag = 1
        P%PureEFTmodelGamma3 = 2
        P%PureEFTHorndeski   = .true.
        EFTparam => P%EFTGamma30
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Gamma3_Hor_linear.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done EFT Gamma3 Horndeski linear in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        write(*,*) 'Pure EFT Horndeski one parameter stability check done.'
        write(*,'(a,F9.2,a,F9.5,a)') ' Total time: ', total_time, &
                                       &' sec. Average time per call: ', total_time/total_number, ' sec.'

    end subroutine SamplePureEFT_Hordenski_1D

    ! ----------------------------

    subroutine Sample_EFT_Horava_1D(NumPoints, IntMin, IntMax, OutRoot, MaxRecursions)
        ! Subroutine to test at once every 1 parameter Horava extensions of the standard model.
        implicit none

        integer , intent(in) :: NumPoints
        real(dl), intent(in) :: IntMin, IntMax

        integer, optional :: MaxRecursions

        character(LEN=100), intent(in) :: OutRoot
        character(LEN=100)             :: OutDir
        real(dl), allocatable :: res(:,:)

        ! timings:
        real(dl) :: t1, t2, total_time = 0._dl
        integer  :: total_number = 0

        ! 1) Set MaxRecursions
        if (.not.present(MaxRecursions)) MaxRecursions = 5
        ! 2) Set the output directory
        OutDir = TRIM(OutRoot)//'11_Horava_1D/'
        call system('mkdir '//OutDir)
        ! 3) Some feedback
        write(*,*) 'Starting Horava one parameter stability check:'
        write(*,*) ' results files placed in: '//TRIM(OutDir)

        ! 1) Horava constant xi
        call EFT_Clean_Params(P)
        P%EFTflag = 4
        P%FullMappingEFTmodel = 1
        EFTparam => P%Horava_xi
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Horava_xi.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Horava xi in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 2) Horava constant lambda
        call EFT_Clean_Params(P)
        P%EFTflag = 4
        P%FullMappingEFTmodel = 1
        EFTparam => P%Horava_lambda
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'2_Horava_lambda.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Horava lambda in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 3) Horava constant eta
        call EFT_Clean_Params(P)
        P%EFTflag = 4
        P%FullMappingEFTmodel = 1
        EFTparam => P%Horava_eta
        t1 = omp_get_wtime()
        call Adaptive_Sampling1D(Stability1D, IntMin, IntMax, NumPoints, res, Recursions = MaxRecursions)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'3_Horava_eta.dat',100)
        write (100,'(E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Horava eta in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        write(*,*) 'Horava one parameter stability check done.'
        write(*,'(a,F9.2,a,F9.5,a)') ' Total time: ', total_time, &
                                       &' sec. Average time per call: ', total_time/total_number, ' sec.'

    end subroutine Sample_EFT_Horava_1D

    ! ----------------------------

    subroutine Sample_EFT_Horava_2D(NumPoints, IntMin1, IntMax1, IntMin2, IntMax2, OutRoot)
        ! Subroutine to test at once every 2 parameter Horava extension of the standard model.
        implicit none

        integer , intent(in) :: NumPoints
        real(dl), intent(in) :: IntMin1, IntMax1, IntMin2, IntMax2
        character(LEN=100), intent(in) :: OutRoot
        character(LEN=100)             :: OutDir

        real(dl), allocatable :: res(:,:)

        ! timings:
        real(dl) :: t1, t2, total_time = 0._dl
        integer  :: total_number = 0

        ! 1) Set the output directory
        OutDir = TRIM(OutRoot)//'12_Horava_2D/'
        call system('mkdir '//OutDir)
        ! 2) Some feedback
        write(*,*) 'Starting Horava two parameter:'
        write(*,*) ' results files placed in: '//TRIM(OutDir)

        ! 1) Horava xi and lambda:
        call EFT_Clean_Params(P)
        P%EFTflag = 4
        P%FullMappingEFTmodel = 1
        EFTparam  => P%Horava_xi
        EFTparam2 => P%Horava_lambda
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'1_Horava_xi_lambda.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Horava xi-lambda in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 2) Horava lambda and eta:
        call EFT_Clean_Params(P)
        P%EFTflag = 4
        P%FullMappingEFTmodel = 1
        EFTparam  => P%Horava_lambda
        EFTparam2 => P%Horava_eta
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'2_Horava_lambda_eta.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Horava lambda-eta in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        ! 3) Horava xi and eta:
        call EFT_Clean_Params(P)
        P%EFTflag = 4
        P%FullMappingEFTmodel = 1
        EFTparam  => P%Horava_xi
        EFTparam2 => P%Horava_eta
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'3_Horava_xi_eta.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Horava xi-eta in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        ! 4) Horava xi and lambda with Solar System Free: IW
        call EFT_Clean_Params(P)
        P%EFTflag = 4
        P%FullMappingEFTmodel = 1
        P%HoravaSolarSystem = .true.
        EFTparam  => P%Horava_xi
        EFTparam2 => P%Horava_lambda
        t1 = omp_get_wtime()
        call Grid_Sampling2D(Stability2D, IntMin1, IntMax1, IntMin2, IntMax2, NumPoints, res)
        t2 = omp_get_wtime()
        call CreateTxtFile(TRIM(OutDir)//'4_Horava_xi_lambda_SolSyst.dat',100)
        write (100,'(E15.5,E15.5,E15.5)') res
        close(100)
        write(*,'(a,F9.2,a,E10.2E2,a)') ' Done Horava xi-lambda SolarSystem in: ', t2-t1, &
            &' sec. Timing per call: ', (t2-t1)/(SIZE(res)/2), ' sec.'

        total_time = total_time +t2-t1
        total_number = total_number + SIZE(res)/2

        write(*,*) 'Horava 2D check done.'
        write(*,'(a,F9.2,a,F9.5,a)') ' Total time: ', total_time, &
                                       &' sec. Average time per call: ', total_time/total_number, ' sec.'

    end subroutine Sample_EFT_Horava_2D

end module Systematic_Model_Sampling

! -------------------------------------------------------------------------------------------------

program EFTStabilityInPSpace
    use CAMB
    use EFTinitialization
    use StabilityFunctions
    use SamplersSubroutines
    use Systematic_Model_Sampling
    implicit none

    real(dl), allocatable :: res(:,:)
    integer :: i
    integer :: x

    real(dl) :: x1,x2, y1, y2
    integer  :: counter, MaxRecursion, num, rec

    character(LEN=100) OutRoot

    ! Gets results root from command line

    OutRoot = ''
    if (GetParamCount() /= 0) then
        OutRoot = GetParam(1)
    else
        OutRoot = 'Stability_Results'
    end if

    OutRoot = TRIM(OutRoot)//'/'

    write(*,*) 'Results files placed in: '//TRIM(OutRoot)
    ! create the results directory:
    call system('mkdir '//(OutRoot))

    ! 0) Set some useful stuff:
    feedbackLevel = 0

    ! 1) Set fixed cosmological parameters:
    call CAMB_SetDefParams(P)

    P%tcmb    = 2.7255
    P%omegab  = 0.046
    P%omegac  = 0.228
    P%omegan  = 0.0
    P%H0      = 70
    P%omegav  = 1-(P%omegab+P%omegac+P%omegan)

    ! stability flags:
    P%EFT_mathematical_stability = .false.
    P%EFT_physical_stability     = .true.
    P%EFTAdditionalPriors        = .true.
    P%MinkowskyPriors            = .true.

    ! 2) fix other cosmological parameters that are not used unless testCls is set to true:

    ! This will call CAMB_GetResults for all the parameters that pass the stability check. Usefull for testing.
    testCls = .false.

    if (testCls) then
        ! This will print sigma 8 for all the models that pass the stability check and that call Cls.
        ! Usefull to spot strange behaviour that do not result in crashes.
        showS8  = .true.
        ! other parameters
        EFTCAMBuseinCOSMOMC = 0

        P%yhe    = 0.24
        P%InitPower%nn     = 1 !number of initial power spectra
        P%InitPower%an(1)  = 1 !scalar spectral index
        P%InitPower%ant(1) = 0 !tensor spectra index
        P%OutputNormalization = outNone

        P%WantScalars = .true.
        P%WantTensors = .true.
        P%DoLensing   = .false.
        P%WantTransfer= .true.
        P%NonLinear   = .false.

        P%Max_l             = 3000
        P%Max_l_tensor      = 1500
        P%Max_eta_k_tensor  = 3000
        P%AccuratePolarization = .true.

        P%Transfer%high_precision=.false.
        P%Transfer%kmax=0.5
        P%Transfer%k_per_logint=3
        P%Transfer%num_redshifts=1
        P%Transfer%redshifts(1)=0
        P%Transfer%PK_num_redshifts=1
        P%Transfer%PK_redshifts(1)=0
        call Transfer_SortAndIndexRedshifts(P%Transfer)

    end if

    ! 3) Do the test of different models:

    if (testCls) then
        num = 10
        rec = 2
    else
        num = 50
        rec = 5
    end if

    call SamplePureEFT1D(num, -2._dl, 2._dl, OutRoot, MaxRecursions = rec)
    call SamplePureEFTwCDM_2D(num, -2._dl, 2._dl, -2._dl, 0._dl, OutRoot)
    call SamplePureEFT_PowerLaw(num, -2._dl, 2._dl, 0._dl, 6._dl, OutRoot)
    call SamplePureEFT_Exponential(num, -2._dl, 2._dl, 0._dl, 6._dl, OutRoot)
    call SampleDesigner1D(num, OutRoot, MaxRecursions = rec)
    call SampleDesigner2D(num, OutRoot)
    call SampleRPH1D(num, -2._dl, 2._dl, OutRoot, MaxRecursions = rec)
    call SampleRPH2D(num, -2._dl, 2._dl, 0._dl, 6._dl, OutRoot)
    call SampleRPH2D_constant(num, -2._dl, 2._dl, -2._dl, 2._dl, OutRoot)
    call SamplePureEFT_Hordenski_1D(num, -2._dl, 2._dl, OutRoot, MaxRecursions = rec)
    call Sample_EFT_Horava_1D(num, -2._dl, 2._dl, OutRoot, MaxRecursions = rec)
    call Sample_EFT_Horava_2D(num, -2._dl, 2._dl, -2._dl, 2._dl, OutRoot)

end program EFTStabilityInPSpace

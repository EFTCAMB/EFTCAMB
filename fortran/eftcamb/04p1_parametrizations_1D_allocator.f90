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

!> @file 04p1_parametrizations_1D_allocator.f90
!! This file contains an interface layer that imports all modules with
!! implementation of abstract parametrization. It also has a function to
!! allocate a generic parametrized function given a model selection flag.

!----------------------------------------------------------------------------------------
!> This module contains an interface layer that imports all modules with
!! implementation of abstract parametrization. It also has a function to
!! allocate a generic parametrized function given a model selection flag.

!> @author Marco Raveri

module EFTCAMB_parametrizations_1D

    ! general modules:
    use EFT_def
    use EFTCAMB_abstract_parametrizations_1D

    ! specific implementations:
    use EFTCAMB_neutral_parametrization_1D
    use EFTCAMB_constant_parametrization_1D
    use EFTCAMB_linear_parametrizations_1D
    use EFTCAMB_power_law_parametrizations_1D
    use EFTCAMB_exponential_parametrizations_1D
    use EFTCAMB_CPL_parametrizations_1D
    use EFTCAMB_JBP_parametrizations_1D
    use EFTCAMB_turning_point_parametrizations_1D
    use EFTCAMB_taylor_parametrizations_1D
    use EFTCAMB_taylor_series_1D
    use EFTCAMB_pade_series_1D
    use EFTCAMB_fourier_parametrizations_1D
    use EFTCAMB_exponential_parametrizations_2_1D
    use EFTCAMB_double_exponential_parametrizations_1D
    use EFTCAMB_cosine_parametrizations_1D
    use EFTCAMB_axion_parametrizations_1D
    use EFTCAMB_step_parametrizations_1D
    use EFTCAMB_steplog_parametrizations_1D
    use EFTCAMB_spline_parametrizations_1D

    implicit none

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates a parametrized 1d function.
    subroutine allocate_parametrized_1D_function( in_function, model_flag, function_name, function_label, eft_error, feedback )

        implicit none

        class( parametrized_function_1D ), allocatable :: in_function    !< The input function to be allocated.
        integer, intent(in)                            :: model_flag     !< the model flag to use to decide allocation
        character(*)                                   :: function_name  !< name of the function
        character(*)                                   :: function_label !< label of the function
        integer                                        :: eft_error      !< error code: 0 all fine, 1 initialization failed
        integer                                        :: feedback       !< this function will print feedback to screen if greater than zero

        ! check allocation and deallocate in case:
        if ( allocated(in_function) ) deallocate(in_function)

        ! perform allocation:
        select case ( model_flag )
            case(0)
                allocate( zero_parametrization_1D::in_function )
            case(1)
                allocate( constant_parametrization_1D::in_function )
            case(2)
                allocate( linear_parametrization_1D::in_function )
            case(3)
                allocate( power_law_parametrization_1D::in_function )
                call in_function%set_param_names( [adjustr(adjustl(function_name))//'0  ', adjustr(adjustl(function_name))//'Exp'], &
                                                & [adjustr(adjustl(function_label))//'_0', adjustr(adjustl(function_label))//'_n'] )
            case(4)
                allocate( exponential_parametrization_1D::in_function )
                call in_function%set_param_names( [adjustr(adjustl(function_name))//'0  ', adjustr(adjustl(function_name))//'Exp'], &
                                                & [adjustr(adjustl(function_label))//'_0', adjustr(adjustl(function_label))//'_n'] )
            case(5)
                allocate( taylorseries_parametrization_1D::in_function )
            case(6)
                allocate( padeseries_parametrization_1D::in_function )
            case(7)
                allocate( fourier_parametrization_1D::in_function )
            case(8)
                allocate( steplog_parametrization_1D::in_function )
            case(9)
                allocate( spline_parametrization_1D::in_function )
            case default
                if (feedback > 0) then
                    write(*,'(a,I3)') 'No model corresponding to flag =', model_flag
                    write(*,'(a,a)')  'For function =', function_name
                    write(*,'(a)')    'Please select an appropriate model.'
                end if
                eft_error = 1
                return
        end select

        ! set names:
        call in_function%set_name ( function_name, function_label )

    end subroutine allocate_parametrized_1D_function

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_parametrizations_1D

!----------------------------------------------------------------------------------------

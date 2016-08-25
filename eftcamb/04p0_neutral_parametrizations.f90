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

!> @file 04p0_neutral_parametrizations.f90
!! This file contains the definition of neutral parametrizations that can be used when
!! parametrized functions should take their GR value.


!----------------------------------------------------------------------------------------
!> This module contains the definition of neutral parametrizations that can be used when
!! parametrized functions should take their GR value.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_neutral_parametrization_1D

    use precision
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the constant function parametrization that returns all zero values.
    !! Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: zero_parametrization_1D

    contains

        ! initialization:
        procedure :: init                  => ZeroParametrized1DInitialize          !< subroutine that initializes the zero parametrization.

        ! evaluation procedures:
        procedure :: value                 => ZeroParametrized1DValue               !< function that returns the value of the zero function.
        procedure :: first_derivative      => ZeroParametrized1DFirstDerivative     !< function that returns the first derivative of the zero function.
        procedure :: second_derivative     => ZeroParametrized1DSecondDerivative    !< function that returns the second derivative of the zero function.
        procedure :: third_derivative      => ZeroParametrized1DThirdDerivative     !< function that returns the third derivative of the zero function.
        procedure :: integral              => ZeroParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type zero_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the constant function parametrization that returns the correct values
    !! to be used as w_DE for a LCDM expansion history. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: wDE_LCDM_parametrization_1D

    contains

        ! initialization:
        procedure :: init                  => wDELCDMParametrized1DInitialize          !< subroutine that initializes the LCDM wDE parametrization.

        ! evaluation procedures:
        procedure :: value                 => wDELCDMParametrized1DValue               !< function that returns the value of the LCDM wDE function.
        procedure :: first_derivative      => wDELCDMParametrized1DFirstDerivative     !< function that returns the first derivative of the LCDM wDE function.
        procedure :: second_derivative     => wDELCDMParametrized1DSecondDerivative    !< function that returns the second derivative of the LCDM wDE function.
        procedure :: third_derivative      => wDELCDMParametrized1DThirdDerivative     !< function that returns the third derivative of the LCDM wDE function.
        procedure :: integral              => wDELCDMParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type wDE_LCDM_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the constant zero function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the zero parametrization.
    subroutine ZeroParametrized1DInitialize( self, name, latexname )

        implicit none

        class(zero_parametrization_1D)  :: self       !< the base class
        character(*), intent(in)        :: name       !< the name of the function
        character(*), intent(in)        :: latexname  !< the latex name of the function

        ! store the name of the function:
        self%name             = TRIM( name )
        ! store the latex name of the function:
        self%name_latex       = TRIM( latexname )
        ! initialize the number of parameters:
        self%parameter_number = 0

    end subroutine ZeroParametrized1DInitialize

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the zero function.
    function ZeroParametrized1DValue( self, x, eft_cache )

        implicit none

        class(zero_parametrization_1D)                     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ZeroParametrized1DValue                             !< the output value

        ZeroParametrized1DValue = 0._dl

    end function ZeroParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the first derivative of the zero function.
    function ZeroParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(zero_parametrization_1D)                     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ZeroParametrized1DFirstDerivative                   !< the output value

        ZeroParametrized1DFirstDerivative = 0._dl

    end function ZeroParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the zero function.
    function ZeroParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(zero_parametrization_1D)                     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ZeroParametrized1DSecondDerivative                  !< the output value

        ZeroParametrized1DSecondDerivative = 0._dl

    end function ZeroParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the zero function.
    function ZeroParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(zero_parametrization_1D)                     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ZeroParametrized1DThirdDerivative                   !< the output value

        ZeroParametrized1DThirdDerivative = 0._dl

    end function ZeroParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the strange integral that we need for w_DE.
    function ZeroParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(zero_parametrization_1D)                     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ZeroParametrized1DIntegral                          !< the output value

        ZeroParametrized1DIntegral = 1._dl/x

    end function ZeroParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the LCDM w_DE.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the LCDM wDE parametrization.
    subroutine wDELCDMParametrized1DInitialize( self, name, latexname )

        implicit none

        class(wDE_LCDM_parametrization_1D)  :: self       !< the base class
        character(*), intent(in)            :: name       !< the name of the function
        character(*), intent(in)            :: latexname  !< the latex name of the function

        ! store the name of the function:
        self%name             = TRIM( name )
        ! store the latex name of the function:
        self%name_latex       = TRIM( latexname )
        ! initialize the number of parameters:
        self%parameter_number = 0

    end subroutine wDELCDMParametrized1DInitialize

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the LCDM wDE function.
    function wDELCDMParametrized1DValue( self, x, eft_cache )

        implicit none

        class(wDE_LCDM_parametrization_1D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: wDELCDMParametrized1DValue                          !< the output value

        wDELCDMParametrized1DValue = -1._dl

    end function wDELCDMParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the first derivative of the LCDM wDE function.
    function wDELCDMParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(wDE_LCDM_parametrization_1D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: wDELCDMParametrized1DFirstDerivative                !< the output value

        wDELCDMParametrized1DFirstDerivative = 0._dl

    end function wDELCDMParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the LCDM wDE function.
    function wDELCDMParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(wDE_LCDM_parametrization_1D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: wDELCDMParametrized1DSecondDerivative               !< the output value

        wDELCDMParametrized1DSecondDerivative = 0._dl

    end function wDELCDMParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the LCDM wDE function.
    function wDELCDMParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(wDE_LCDM_parametrization_1D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: wDELCDMParametrized1DThirdDerivative                !< the output value

        wDELCDMParametrized1DThirdDerivative = 0._dl

    end function wDELCDMParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the strange integral that we need for w_DE.
    function wDELCDMParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(wDE_LCDM_parametrization_1D)                 :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: wDELCDMParametrized1DIntegral                       !< the output value

        wDELCDMParametrized1DIntegral = x**2

    end function wDELCDMParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_neutral_parametrization_1D

!----------------------------------------------------------------------------------------

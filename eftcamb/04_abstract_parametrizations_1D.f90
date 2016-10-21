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

!> @file 04_abstract_parametrizations_1D.f90
!! This file contains the abstract class for generic parametrizations for 1D functions
!! that are used by several models in EFTCAMB. When there is a free function
!! in EFT it should be declared as a class inheriting from parametrized_function_1D.
!! This guarantees maximum performances as well as maximum flexibility.


!----------------------------------------------------------------------------------------
!> This module contains the abstract class for generic parametrizations for 1D functions
!! that are used by several models in EFTCAMB. When there is a free function
!! in EFT it should be declared as a class inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_abstract_parametrizations_1D

    use precision
    use AMLutils
    use IniFile
    use EFT_def
    use EFTCAMB_mixed_algorithms
    use EFTCAMB_cache

    implicit none

    private

    public parametrized_function_1D

    !----------------------------------------------------------------------------------------
    !> This is the abstract type for parametrized functions. As a rule, when there is a
    !! free function in EFT it should be declared as a class inheriting from parametrized_function_1D.
    !! This guarantees maximum performances as well as maximum flexibility.
    type, abstract :: parametrized_function_1D

        integer                                     :: parameter_number  !< number of parameters defining the parametrized function.
        character(len=:), allocatable               :: name              !< name of the function.
        character(len=:), allocatable               :: name_latex        !< latex name of the function.
        type(string)    , allocatable, dimension(:) :: param_names       !< array of strings with the names of the parameters.
        type(string)    , allocatable, dimension(:) :: param_names_latex !< array of strings with the latex version of the parameter names.

    contains

        ! initialization procedures:
        procedure( ParametrizedFunction1DSetParamNumber   ), deferred :: set_param_number  !< subroutine that sets the number of parameters of the parametrized function.
        procedure :: param_number          => ParametrizedFunction1DParamNumber            !< function that returns the number of parameters of the parametrized function.
        procedure :: set_name              => ParametrizedFunction1DSetName                !< subroutine that sets the name of the parametrized function.
        procedure :: set_param_names       => ParametrizedFunction1DSetParamNames          !< subroutine that sets the parameter names starting from two array of characters, containing the names and the latex names of the function parameters.
        procedure :: init_from_file        => ParametrizedFunction1DInitFromFile           !< subroutine that reads a Ini file looking for the parameters of the function.
        procedure( ParametrizedFunction1DInitParams       ), deferred :: init_parameters   !< subroutine that initializes the function parameters based on the values found in an input array.
        ! utility functions:
        procedure :: feedback              => ParametrizedFunction1DFeedback               !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => ParametrizedFunction1DParameterNames         !< subroutine that returns the i-th parameter name of the function.
        procedure :: parameter_names_latex => ParametrizedFunction1DParameterNamesLatex    !< subroutine that returns the i-th parameter name of the function in latex format.
        procedure( ParametrizedFunction1DParameterValues  ), deferred :: parameter_value   !< subroutine that returns the value of the function i-th parameter.
        ! evaluation procedures:
        procedure( ParametrizedFunction1DValue            ), deferred :: value             !< function that returns the value of the function. The EFTCAMB cache is passed as an optional argument in case the parametrization uses some background quantity.
        procedure( ParametrizedFunction1DFirstDerivative  ), deferred :: first_derivative  !< function that returns the first derivative of the function. The EFTCAMB cache is passed as an optional argument in case the parametrization uses some background quantity.
        procedure( ParametrizedFunction1DSecondDerivative ), deferred :: second_derivative !< function that returns the second derivative of the function. The EFTCAMB cache is passed as an optional argument in case the parametrization uses some background quantity.
        procedure( ParametrizedFunction1DThirdDerivative  ), deferred :: third_derivative  !< function that returns the third derivative of the function. The EFTCAMB cache is passed as an optional argument in case the parametrization uses some background quantity.
        procedure( ParametrizedFunction1DIntegral         ), deferred :: integral          !< function that returns the strange integral that we need for w_DE. The EFTCAMB cache is passed as an optional argument in case the parametrization uses some background quantity.

    end type parametrized_function_1D

    ! ---------------------------------------------------------------------------------------------
    ! parametrized_function_1D abstract interfaces: these are all the function procedures
    ! that the user HAS to override when writing its own parametrized 1D function.
    ! ---------------------------------------------------------------------------------------------

    abstract interface

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that sets the number of parameters of the parametrized function.
        subroutine ParametrizedFunction1DSetParamNumber( self )
            use precision
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D)  :: self       !< the base class
        end subroutine ParametrizedFunction1DSetParamNumber

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that reads a Ini file looking for the parameters of the function.
        subroutine ParametrizedFunction1DInitParams( self, array )
            use precision
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D)                        :: self   !< the base class.
            real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.
        end subroutine ParametrizedFunction1DInitParams

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that returns the value of the function i-th parameter.
        subroutine ParametrizedFunction1DParameterValues( self, i, value )
            use precision
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D) :: self   !< the base class.
            integer , intent(in)            :: i      !< input number of the parameter
            real(dl), intent(out)           :: value  !< output value of the parameter
        end subroutine ParametrizedFunction1DParameterValues

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the function. The EFTCAMB cache is passed as an optional
        !! argument in case the parametrization uses some background quantity.
        function ParametrizedFunction1DValue( self, x, eft_cache )
            use    precision
            use    EFTCAMB_cache
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D)                    :: self      !< the base class
            real(dl), intent(in)                               :: x         !< the input scale factor
            type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
            real(dl) :: ParametrizedFunction1DValue                         !< the output value
        end function ParametrizedFunction1DValue

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the first derivative of the function
        !! with respect to the scale factor. The EFTCAMB cache is passed as an optional
        !! argument in case the parametrization uses some background quantity.
        function ParametrizedFunction1DFirstDerivative( self, x, eft_cache )
            use    precision
            use    EFTCAMB_cache
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D)                    :: self      !< the base class
            real(dl), intent(in)                               :: x         !< the input scale factor
            type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
            real(dl) :: ParametrizedFunction1DFirstDerivative               !< the output value
        end function ParametrizedFunction1DFirstDerivative

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the second derivative of the function
        !! with respect to the scale factor. The EFTCAMB cache is passed as an optional
        !! argument in case the parametrization uses some background quantity.
        function ParametrizedFunction1DSecondDerivative( self, x, eft_cache )
            use    precision
            use    EFTCAMB_cache
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D)                    :: self      !< the base class
            real(dl), intent(in)                               :: x         !< the input scale factor
            type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
            real(dl) :: ParametrizedFunction1DSecondDerivative              !< the output value
        end function ParametrizedFunction1DSecondDerivative

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the value of the third derivative of the function
        !! with respect to the scale factor. The EFTCAMB cache is passed as an optional
        !! argument in case the parametrization uses some background quantity.
        function ParametrizedFunction1DThirdDerivative( self, x, eft_cache )
            use    precision
            use    EFTCAMB_cache
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D)                    :: self      !< the base class
            real(dl), intent(in)                               :: x         !< the input scale factor
            type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
            real(dl) :: ParametrizedFunction1DThirdDerivative               !< the output value
        end function ParametrizedFunction1DThirdDerivative

        ! ---------------------------------------------------------------------------------------------
        !> Function that returns the integral of the function, as defined in the notes.
        !! The EFTCAMB cache is passed as an optional argument in case the parametrization uses some background quantity.
        function ParametrizedFunction1DIntegral( self, x, eft_cache )
            use    precision
            use    EFTCAMB_cache
            import parametrized_function_1D
            implicit none
            class(parametrized_function_1D)                    :: self      !< the base class
            real(dl), intent(in)                               :: x         !< the input scale factor
            type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
            real(dl) :: ParametrizedFunction1DIntegral                      !< the output value
        end function ParametrizedFunction1DIntegral

        ! ---------------------------------------------------------------------------------------------

    end interface

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the number of parameters of the parametrized function.
    function ParametrizedFunction1DParamNumber( self )

        implicit none

        class(parametrized_function_1D)    :: self       !< the base class
        integer  :: ParametrizedFunction1DParamNumber    !< the number of parameters of the parametrized function

        ParametrizedFunction1DParamNumber = self%parameter_number

    end function ParametrizedFunction1DParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the name of the parametrized function.
    subroutine ParametrizedFunction1DSetName( self, name, latexname )

        implicit none

        class(parametrized_function_1D)    :: self       !< the base class
        character(*), intent(in)           :: name       !< the name of the function
        character(*), intent(in), optional :: latexname  !< the optional latex name of the function

        ! ensure that the number of parameters is properly associated:
        call self%set_param_number()

        ! store the name of the function:
        self%name             = TRIM( name      )
        ! store the latex name of the function:
        if ( present(latexname) ) then
            self%name_latex   = TRIM( latexname )
        else
            self%name_latex   = TRIM( name      )
        end if

    end subroutine ParametrizedFunction1DSetName

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the parameter names starting from two array of characters,
    !! containing the names and the latex names of the function parameters.
    subroutine ParametrizedFunction1DSetParamNames( self, param_names, param_names_latex )

        implicit none

        class(parametrized_function_1D)                   :: self              !< the base class
        character(*), intent(in), dimension(:)            :: param_names       !< an array of strings containing the names of the function parameters
        character(*), intent(in), dimension(:), optional  :: param_names_latex !< an array of strings containing the latex names of the function parameters

        character(len=:), allocatable :: string_array(:)

        integer  :: num_params, ind

        ! ensure that the number of parameters is properly associated:
        call self%set_param_number()

        ! check the number of parameters:
        num_params = self%param_number()
        if (  num_params /= size(param_names) ) then
            write(*,*) 'In parametrized_function_1D:', self%name
            write(*,*) 'Length of param_names and number of parameters do not coincide.'
            write(*,*) 'Parameter number:', num_params
            write(*,*) 'Size of the param_names array:', size(param_names)
            call MpiStop('EFTCAMB error')
        end if
        if ( present(param_names_latex) ) then
            ! check length:
            if (  num_params /= size(param_names_latex) ) then
                write(*,*) 'In parametrized_function_1D:', self%name
                write(*,*) 'Length of param_names_latex and number of parameters do not coincide.'
                write(*,*) 'Parameter number:', self%parameter_number
                write(*,*) 'Size of the param_names array:', size(param_names_latex)
                call MpiStop('EFTCAMB error')
            end if
        end if

        ! allocate self%param_names and self%param_names_latex:
        if ( allocated(self%param_names) ) deallocate(self%param_names)
        allocate( self%param_names(num_params) )
        if ( allocated(self%param_names_latex) ) deallocate(self%param_names_latex)
        allocate( self%param_names_latex(num_params) )

        ! store the parameter names and latex param names:
        do ind=1, num_params
            self%param_names(ind)%string  = param_names(ind)
            if ( present(param_names_latex) ) then
                self%param_names_latex(ind)%string = param_names_latex(ind)
            else
                self%param_names_latex(ind)%string = param_names(ind)
            end if
        end do

    end subroutine ParametrizedFunction1DSetParamNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads a Ini file looking for the parameters of the function.
    subroutine ParametrizedFunction1DInitFromFile( self, Ini )

        implicit none

        class(parametrized_function_1D)    :: self   !< the base class
        type(TIniFile)                     :: Ini    !< Input ini file

        character(len=EFT_names_max_length)          :: param_name
        real(dl), dimension( self%parameter_number ) :: parameters

        integer  :: i

        ! ensure that the number of parameters is properly associated:
        call self%set_param_number()

        ! read the parameters and store them in a vector:
        do i=1, self%parameter_number
            call self%parameter_names( i, param_name )
            parameters(i) = Ini_Read_Double_File( Ini, TRIM(param_name), 0._dl )
        end do
        ! initialize the function parameters from the vector:
        call self%init_parameters( parameters )

    end subroutine ParametrizedFunction1DInitFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine ParametrizedFunction1DFeedback( self )

        implicit none

        class(parametrized_function_1D)  :: self   !< the base class

        integer                             :: i
        real(dl)                            :: param_value
        character(len=EFT_names_max_length) :: param_name

        if ( self%parameter_number>0 ) then
            write(*,*)     'Parametrized function 1D: ', self%name
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine ParametrizedFunction1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name.
    subroutine ParametrizedFunction1DParameterNames( self, i, name )

        implicit none

        class(parametrized_function_1D) :: self   !< the base class
        integer     , intent(in)        :: i      !< the index of the parameter
        character(*), intent(out)       :: name   !< the output name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'In parametrized_function_1D:', self%name
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the parameter name:
        if ( allocated(self%param_names) ) then
            name = self%param_names(i)%string
        else
            name = TRIM(self%name)//integer_to_string(i-1)
        end if

    end subroutine ParametrizedFunction1DParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name.
    subroutine ParametrizedFunction1DParameterNamesLatex( self, i, latexname )

        implicit none

        class(parametrized_function_1D) :: self       !< the base class
        integer     , intent(in)        :: i          !< the index of the parameter
        character(*), intent(out)       :: latexname  !< the output latex name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'In parametrized_function_1D:', self%name
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the parameter name:
        if ( allocated(self%param_names_latex) ) then
            latexname = self%param_names_latex(i)%string
        else
            latexname = TRIM(self%name)//'_'//integer_to_string(i-1)
        end if

    end subroutine ParametrizedFunction1DParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_abstract_parametrizations_1D

!----------------------------------------------------------------------------------------

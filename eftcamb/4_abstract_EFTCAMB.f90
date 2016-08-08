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

!> @file 4_abstract_EFTCAMB.f90
!! This file contains the abstract definition of all the places where EFTCAMB interacts
!! with CAMB. All EFTCAMB models should inherit from this class.


!----------------------------------------------------------------------------------------
!> This module contains the abstract definition of all the places where EFTCAMB interacts
!! with CAMB. All EFTCAMB models should inherit from this class.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_abstract_model

    use precision
    use IniFile

    implicit none

    !----------------------------------------------------------------------------------------
    !> This is the abstract type for EFTCAMB models. As a rule, when there is a
    !! new model it should be declared as a class inheriting from EFTCAMB_model.
    !! This guarantees maximum performances as well as maximum flexibility.
    type EFTCAMB_model

        integer                       :: parameter_number !< number of parameters defining the model
        character(len=:), allocatable :: name             !< name of the model
        character(len=:), allocatable :: name_latex       !< latex name of the model

    contains

        ! utility functions:
        procedure :: model_name            => EFTCAMBModelName                !< subroutine that returns the name of the EFTCAMB model
        procedure :: parameter_names       => EFTCAMBModelParameterNames      !< subroutine that returns the i-th parameter name of the model
        procedure :: parameter_names_latex => EFTCAMBModelParameterNamesLatex !< subroutine that returns the i-th parameter name of the model
        procedure :: parameter_values      => EFTCAMBModelParameterValues     !< subroutine that returns the i-th parameter value
        ! initialization of the model:
        procedure :: read_parameters_file  => EFTCAMBModelReadParametersFromFile !< subroutine that reads the parameters of the model from file

    end type EFTCAMB_model

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the name of the EFTCAMB model
    subroutine EFTCAMBModelName( self, name )

        implicit none

        class(EFTCAMB_model)       :: self   !< the base class
        character(*), intent(out)  :: name   !< the output name of the model

        ! check that the model name is initialized:
        if ( .not. allocated(self%name) ) then
            write(*,*) 'ERROR, the model name is not set. Most likely you do not have run properly model initialization.'
            stop
        end if
        ! return the name:
        name = self%name

    end subroutine EFTCAMBModelName

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBModelParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_model)      :: self   !< the base class
        integer     , intent(in)  :: i      !< the index of the parameter
        character(*), intent(out) :: name   !< the output name of the i-th parameter

    end subroutine EFTCAMBModelParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBModelParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_model)      :: self       !< the base class
        integer     , intent(in)  :: i          !< The index of the parameter
        character(*), intent(out) :: latexname  !< the output latex name of the i-th parameter

    end subroutine EFTCAMBModelParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBModelParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_model)  :: self   !< the base class
        integer , intent(in)  :: i      !< The index of the parameter
        real(dl), intent(out) :: value  !< the output value of the i-th parameter

    end subroutine EFTCAMBModelParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBModelReadParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_model)  :: self   !< the base class
        type(TIniFile)        :: Ini    !< Input ini file

    end subroutine EFTCAMBModelReadParametersFromFile

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_abstract_model

!----------------------------------------------------------------------------------------

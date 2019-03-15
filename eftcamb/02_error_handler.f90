!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2019 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 02_error_handler.f90
!! This file contains the EFTCAMB error handler that gives traceback informations.

!----------------------------------------------------------------------------------------
!> This module contains the EFTCAMB error handler that gives traceback informations.

!> @author Marco Raveri

module EFTCAMB_error_handler

    implicit none

    private

    !----------------------------------------------------------------------------------------
    ! Internal definitions:
    interface
        subroutine abort() bind(C, name="abort")
        end subroutine
    end interface

    !----------------------------------------------------------------------------------------
    ! public definitions:
    public :: EFTCAMB_error

    !----------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that raises EFTCAMB errors
    subroutine EFTCAMB_error( error_message )

        implicit none

        character ( len = *  ), intent(in), optional :: error_message !< input error message

        ! if there is a message present print it:
        if ( present( error_message ) ) then
            write(*,*) 'EFTCAMB error:'
            write(*,*) error_message
        end if
        ! call abort function that prints traceback infos when compiled in debug mode:
        call abort()

    end subroutine EFTCAMB_error

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_error_handler

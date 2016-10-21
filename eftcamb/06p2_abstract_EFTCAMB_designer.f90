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

!> @file 06p2_abstract_EFTCAMB_designer.f90
!! This file contains the abstract definition of all the places where EFTCAMB interacts
!! with CAMB in case of a designer model. All models implementing a model in which
!! the cosmological background is parametrized  according to some choice of the user
!! should inherit from this class.


!----------------------------------------------------------------------------------------
!> This module contains the abstract definition of all the places where EFTCAMB interacts
!! with CAMB in case of a designer model. All models implementing a model in which
!! the cosmological background is parametrized  according to some choice of the user
!! should inherit from this class.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_abstract_model_designer

    use precision
    use IniFile
    use EFTCAMB_abstract_model

    implicit none

    private

    public EFTCAMB_designer_model

    !----------------------------------------------------------------------------------------
    !> This is the abstract type for EFTCAMB models.
    type, extends ( EFTCAMB_model ), abstract :: EFTCAMB_designer_model

    end type EFTCAMB_designer_model

contains

end module EFTCAMB_abstract_model_designer

!----------------------------------------------------------------------------------------

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

!> @file 5_Pure_EFT_std.f90
!! This file contains the definition of the Pure EFT model in which the EFT is described
!! by six functions of time and w_DE. Please refer to the numerical notes for details.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Pure EFT model in which the EFT is described
!! by six functions of time and w_DE. Please refer to the numerical notes for details.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_pure_EFT_std

    use precision
    use IniFile
    use EFTCAMB_abstract_parametrizations
    use EFTCAMB_abstract_model

    implicit none

    !----------------------------------------------------------------------------------------
    !> This is the pure EFT model with the six functions of time and w_DE.
    type, extends ( EFTCAMB_model ) :: EFTCAMB_std_pure_EFT

        ! the pure EFT functions model selection flags:
        integer  :: PureEFTmodelOmega   !< Model selection flag for Pure EFT Omega
        integer  :: EFTwDE              !< Model selection flag for Pure EFT w DE
        integer  :: PureEFTmodelGamma1  !< Model selection flag for Pure EFT Gamma1
        integer  :: PureEFTmodelGamma2  !< Model selection flag for Pure EFT Gamma2
        integer  :: PureEFTmodelGamma3  !< Model selection flag for Pure EFT Gamma3
        integer  :: PureEFTmodelGamma4  !< Model selection flag for Pure EFT Gamma4
        integer  :: PureEFTmodelGamma5  !< Model selection flag for Pure EFT Gamma5
        integer  :: PureEFTmodelGamma6  !< Model selection flag for Pure EFT Gamma6

        ! the pure EFT functions:
        class( parametrized_function ), allocatable :: PureEFTOmega    !< The pure EFT function Omega
        class( parametrized_function ), allocatable :: PureEFTwDE      !< The pure EFT function w_DE
        class( parametrized_function ), allocatable :: PureEFTGamma1   !< The pure EFT function Gamma1
        class( parametrized_function ), allocatable :: PureEFTGamma2   !< The pure EFT function Gamma2
        class( parametrized_function ), allocatable :: PureEFTGamma3   !< The pure EFT function Gamma3
        class( parametrized_function ), allocatable :: PureEFTGamma4   !< The pure EFT function Gamma4
        class( parametrized_function ), allocatable :: PureEFTGamma5   !< The pure EFT function Gamma5
        class( parametrized_function ), allocatable :: PureEFTGamma6   !< The pure EFT function Gamma6

    contains

        ! utility functions:
        procedure :: feedback              => EFTCAMBPureEFTstdFeedback            !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBPureEFTstdParameterNames      !< subroutine that returns the i-th parameter name of the model
        procedure :: parameter_names_latex => EFTCAMBPureEFTstdParameterNamesLatex !< subroutine that returns the i-th parameter name of the model
        procedure :: parameter_values      => EFTCAMBPureEFTstdParameterValues     !< subroutine that returns the i-th parameter value
        ! initialization of the model:
        procedure :: read_model_selection  => EFTCAMBPureEFTstdReadModelSelectionFromFile !< subroutine that reads the parameters of the model from file

    end type EFTCAMB_std_pure_EFT

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBPureEFTstdFeedback( self )

        implicit none

        class(EFTCAMB_std_pure_EFT)  :: self   !< the base class

        write(*,*)
        write(*,'(a)')    '   Model               =', self%name
        ! print model functions informations:
        write(*,*)
        write(*,'(a,I3)') '   PureEFTmodelOmega   =', self%PureEFTmodelOmega
        call self%PureEFTOmega%feedback()
        write(*,'(a,I3)') '   EFTwDE              =', self%EFTwDE
        call self%PureEFTwDE%feedback()
        write(*,'(a,I3)') '   PureEFTmodelGamma1  =', self%PureEFTmodelGamma1
        call self%PureEFTGamma1%feedback()
        write(*,'(a,I3)') '   PureEFTmodelGamma2  =', self%PureEFTmodelGamma2
        call self%PureEFTGamma2%feedback()
        write(*,'(a,I3)') '   PureEFTmodelGamma3  =', self%PureEFTmodelGamma3
        call self%PureEFTGamma3%feedback()
        write(*,'(a,I3)') '   PureEFTmodelGamma4  =', self%PureEFTmodelGamma4
        call self%PureEFTGamma4%feedback()
        write(*,'(a,I3)') '   PureEFTmodelGamma5  =', self%PureEFTmodelGamma5
        call self%PureEFTGamma5%feedback()
        write(*,'(a,I3)') '   PureEFTmodelGamma6  =', self%PureEFTmodelGamma6
        call self%PureEFTGamma6%feedback()

    end subroutine EFTCAMBPureEFTstdFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

    end subroutine EFTCAMBPureEFTstdParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self       !< the base class
        integer     , intent(in)    :: i          !< The index of the parameter
        character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

    end subroutine EFTCAMBPureEFTstdParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBPureEFTstdParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        integer , intent(in)        :: i      !< The index of the parameter
        real(dl), intent(out)       :: value  !< the output value of the i-th parameter

    end subroutine EFTCAMBPureEFTstdParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBPureEFTstdReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_std_pure_EFT) :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        self%PureEFTmodelOmega  = Ini_Read_Int_File( Ini, 'PureEFTmodelOmega'  , 0 )
        self%EFTwDE             = Ini_Read_Int_File( Ini, 'EFTwDE'             , 0 )
        self%PureEFTmodelGamma1 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma1' , 0 )
        self%PureEFTmodelGamma2 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma2' , 0 )
        self%PureEFTmodelGamma3 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma3' , 0 )
        self%PureEFTmodelGamma4 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma4' , 0 )
        self%PureEFTmodelGamma5 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma5' , 0 )
        self%PureEFTmodelGamma6 = Ini_Read_Int_File( Ini, 'PureEFTmodelGamma6' , 0 )

    end subroutine EFTCAMBPureEFTstdReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_pure_EFT_std

!----------------------------------------------------------------------------------------

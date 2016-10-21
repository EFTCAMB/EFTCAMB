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

!> @file 09p2_Designer_mc_quintessence.f90
!! This file contains the definition of the designer minimally coupled quintessence
!! model.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the designer minimally coupled quintessence
!! model.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_designer_mc_quintessence

    use precision
    use IniFile
    use AMLutils
    use EFTCAMB_neutral_parametrization_1D
    use EFTCAMB_constant_parametrization_1D
    use EFTCAMB_linear_parametrizations_1D
    use EFTCAMB_power_law_parametrizations_1D
    use EFTCAMB_exponential_parametrizations_1D
    use EFTCAMB_CPL_parametrizations_1D
    use EFTCAMB_JBP_parametrizations_1D
    use EFTCAMB_turning_point_parametrizations_1D
    use EFTCAMB_taylor_parametrizations_1D
    use EFTCAMB_pure_EFT_std

    implicit none

    private

    public EFTCAMB_des_mc_quint

    !----------------------------------------------------------------------------------------
    !> This is the pure EFT model with the six functions of time and w_DE.
    type, extends ( EFTCAMB_std_pure_EFT ) :: EFTCAMB_des_mc_quint

    contains

        ! initialization of the model:
        procedure :: allocate_model_selection  => EFTCAMBDesMC5eAllocateModelSelection      !< subroutine that allocates the model selection.
        ! utility functions:
        procedure :: feedback                  => EFTCAMBDesMC5eFeedback                    !< subroutine that prints on the screen feedback information about the model.

    end type EFTCAMB_des_mc_quint

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBDesMC5eAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_des_mc_quint)                       :: self              !< the base class

        ! allocate wDE:
        if ( allocated(self%PureEFTwDE) ) deallocate(self%PureEFTwDE)
        select case ( self%EFTwDE )
            case(0)
                allocate( wDE_LCDM_parametrization_1D::self%PureEFTwDE )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTwDE )
            case(2)
                allocate( CPL_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa'], ['w_0', 'w_a'] )
            case(3)
                allocate( JBP_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTwn'], [ 'w_0', 'w_a', 'n  ' ] )
            case(4)
                allocate( turning_point_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0 ', 'EFTwa ', 'EFTwat'], ['w_0', 'w_a', 'a_t'] )
            case(5)
                allocate( taylor_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTw2', 'EFTw3'], ['w_0', 'w_a', 'w_2', 'w_3'] )
            case default
                write(*,'(a,I3)') 'No model corresponding to EFTwDE =', self%EFTwDE
                write(*,'(a)')    'Please select an appropriate model.'
        end select
        ! allocate Omega:
        if ( allocated(self%PureEFTOmega) ) deallocate(self%PureEFTOmega)
        allocate( zero_parametrization_1D::self%PureEFTOmega )
        ! allocate Gamma1:
        if ( allocated(self%PureEFTGamma1) ) deallocate(self%PureEFTGamma1)
        allocate( zero_parametrization_1D::self%PureEFTGamma1 )
        ! allocate Gamma2:
        if ( allocated(self%PureEFTGamma2) ) deallocate(self%PureEFTGamma2)
        allocate( zero_parametrization_1D::self%PureEFTGamma2 )
        ! allocate Gamma3:
        if ( allocated(self%PureEFTGamma3) ) deallocate(self%PureEFTGamma3)
        allocate( zero_parametrization_1D::self%PureEFTGamma3 )
        ! allocate Gamma4:
        if ( allocated(self%PureEFTGamma4) ) deallocate(self%PureEFTGamma4)
        allocate( zero_parametrization_1D::self%PureEFTGamma4 )
        ! allocate Gamma5:
        if ( allocated(self%PureEFTGamma5) ) deallocate(self%PureEFTGamma5)
        allocate( zero_parametrization_1D::self%PureEFTGamma5 )
        ! allocate Gamma6:
        if ( allocated(self%PureEFTGamma6) ) deallocate(self%PureEFTGamma6)
        allocate( zero_parametrization_1D::self%PureEFTGamma6 )

        ! initialize the names:
        call self%PureEFTOmega%set_name ( 'EFTOmega' , '\Omega'       )
        call self%PureEFTwDE%set_name   ( 'EFTw'     , 'w'            )
        call self%PureEFTGamma1%set_name( 'EFTGamma1', '\gamma^{(1)}' )
        call self%PureEFTGamma2%set_name( 'EFTGamma2', '\gamma^{(2)}' )
        call self%PureEFTGamma3%set_name( 'EFTGamma3', '\gamma^{(3)}' )
        call self%PureEFTGamma4%set_name( 'EFTGamma4', '\gamma^{(4)}' )
        call self%PureEFTGamma5%set_name( 'EFTGamma5', '\gamma^{(5)}' )
        call self%PureEFTGamma6%set_name( 'EFTGamma6', '\gamma^{(6)}' )

    end subroutine EFTCAMBDesMC5eAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBDesMC5eFeedback( self )

        implicit none

        class(EFTCAMB_des_mc_quint)  :: self   !< the base class

        ! print general model informations:
        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number

        ! print model functions informations:
        write(*,*)
        if ( self%EFTwDE             /= 0 ) write(*,'(a,I3)') '   EFTwDE              =', self%EFTwDE

        write(*,*)
        ! print functions informations:
        call self%PureEFTwDE%feedback()

    end subroutine EFTCAMBDesMC5eFeedback

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_designer_mc_quintessence

!----------------------------------------------------------------------------------------

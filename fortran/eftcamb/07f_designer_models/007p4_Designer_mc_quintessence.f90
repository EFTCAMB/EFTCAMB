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

!> @file 007p4_Designer_mc_quintessence.f90
!! This file contains the definition of the designer minimally coupled quintessence
!! model.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the designer minimally coupled quintessence
!! model.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_designer_mc_quintessence

    use precision
    use IniObjects
    use EFTCAMB_parametrizations_1D
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
    subroutine EFTCAMBDesMC5eAllocateModelSelection( self, Ini, eft_error )

        implicit none

        class(EFTCAMB_des_mc_quint)   :: self      !< the base class
        type(TIniFile)                :: Ini       !< Input ini file
        integer                       :: eft_error !< error code: 0 all fine, 1 initialization failed

        integer                     :: temp_feedback 

        ! get feedback flag:
        temp_feedback = Ini%Read_Int('feedback_level', 0)

        ! allocate wDE:
        if ( allocated(self%PureEFTwDE) ) deallocate(self%PureEFTwDE)
        select case ( self%EFTwDE )
            case(0)
                allocate( wDE_LCDM_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_name( 'EFTw', 'w' )
            case(1)
                allocate( constant_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_name( 'EFTw', 'w' )
            case(2)
                allocate( CPL_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa'], ['w_0', 'w_a'] )
                call self%PureEFTwDE%set_name( 'EFTw', 'w' )
            case(3)
                allocate( JBP_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTwn'], [ 'w_0', 'w_a', 'n  ' ] )
                call self%PureEFTwDE%set_name( 'EFTw', 'w' )
            case(4)
                allocate( turning_point_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0 ', 'EFTwa ', 'EFTwat'], ['w_0', 'w_a', 'a_t'] )
                call self%PureEFTwDE%set_name( 'EFTw', 'w' )
            case(5)
                allocate( taylor_parametrization_1D::self%PureEFTwDE )
                call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTw2', 'EFTw3'], ['w_0', 'w_a', 'w_2', 'w_3'] )
                call self%PureEFTwDE%set_name( 'EFTw', 'w' )
            case default
                call allocate_parametrized_1D_function( self%PureEFTwDE, self%EFTwDE, 'EFTw', 'w', eft_error, temp_feedback )
                if ( eft_error == 1 ) return
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
        ! allocate Omega_ODE:
        if ( allocated(self%PureEFTOmega_ODE) ) deallocate(self%PureEFTOmega_ODE)
        allocate( zero_parametrization_1D::self%PureEFTOmega_ODE )
        ! allocate Gamma1_ODE:
        if ( allocated(self%PureEFTGamma1_ODE) ) deallocate(self%PureEFTGamma1_ODE)
        allocate( zero_parametrization_1D::self%PureEFTGamma1_ODE )
        ! allocate Gamma2_ODE:
        if ( allocated(self%PureEFTGamma2_ODE) ) deallocate(self%PureEFTGamma2_ODE)
        allocate( zero_parametrization_1D::self%PureEFTGamma2_ODE )
        ! allocate Gamma3_ODE:
        if ( allocated(self%PureEFTGamma3_ODE) ) deallocate(self%PureEFTGamma3_ODE)
        allocate( zero_parametrization_1D::self%PureEFTGamma3_ODE )
        ! allocate Gamma4_ODE:
        if ( allocated(self%PureEFTGamma4_ODE) ) deallocate(self%PureEFTGamma4_ODE)
        allocate( zero_parametrization_1D::self%PureEFTGamma4_ODE )
        ! allocate Gamma5_ODE:
        if ( allocated(self%PureEFTGamma5_ODE) ) deallocate(self%PureEFTGamma5_ODE)
        allocate( zero_parametrization_1D::self%PureEFTGamma5_ODE )
        ! allocate Gamma6_ODE:
        if ( allocated(self%PureEFTGamma6_ODE) ) deallocate(self%PureEFTGamma6_ODE)
        allocate( zero_parametrization_1D::self%PureEFTGamma6_ODE )

        ! initialize the names:
        call self%PureEFTOmega%set_name ( 'EFTOmega' , '\Omega'       )
        call self%PureEFTwDE%set_name   ( 'EFTw'     , 'w'            )
        call self%PureEFTGamma1%set_name( 'EFTGamma1', '\gamma^{(1)}' )
        call self%PureEFTGamma2%set_name( 'EFTGamma2', '\gamma^{(2)}' )
        call self%PureEFTGamma3%set_name( 'EFTGamma3', '\gamma^{(3)}' )
        call self%PureEFTGamma4%set_name( 'EFTGamma4', '\gamma^{(4)}' )
        call self%PureEFTGamma5%set_name( 'EFTGamma5', '\gamma^{(5)}' )
        call self%PureEFTGamma6%set_name( 'EFTGamma6', '\gamma^{(6)}' )
        call self%PureEFTOmega_ODE%set_name ( 'EFTOmega_ODE' , '\Omega_{ODE}'       )
        call self%PureEFTGamma1_ODE%set_name( 'EFTGamma1_ODE', '\gamma^{(1)}_{ODE}' )
        call self%PureEFTGamma2_ODE%set_name( 'EFTGamma2_ODE', '\gamma^{(2)}_{ODE}' )
        call self%PureEFTGamma3_ODE%set_name( 'EFTGamma3_ODE', '\gamma^{(3)}_{ODE}' )
        call self%PureEFTGamma4_ODE%set_name( 'EFTGamma4_ODE', '\gamma^{(4)}_{ODE}' )
        call self%PureEFTGamma5_ODE%set_name( 'EFTGamma5_ODE', '\gamma^{(5)}_{ODE}' )
        call self%PureEFTGamma6_ODE%set_name( 'EFTGamma6_ODE', '\gamma^{(6)}_{ODE}' )

        ! additional initialization of the function:
        call self%PureEFTwDE%init_func_from_file( Ini, eft_error )

    end subroutine EFTCAMBDesMC5eAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBDesMC5eFeedback( self, print_params )

        implicit none

        class(EFTCAMB_des_mc_quint)  :: self         !< the base class
        logical, optional            :: print_params !< optional flag that decised whether to print numerical values
                                                     !! of the parameters.

        ! print general model informations:
        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number

        ! print model functions informations:
        write(*,*)
        if ( self%EFTwDE             /= 0 ) write(*,'(a,I3)') '   EFTwDE              =', self%EFTwDE

        write(*,*)
        ! print functions informations:
        call self%PureEFTwDE%feedback( print_params )

    end subroutine EFTCAMBDesMC5eFeedback

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_designer_mc_quintessence

!----------------------------------------------------------------------------------------

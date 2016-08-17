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

!> @file 05_abstract_EFTCAMB_cache.f90
!! This file contains the definition of a time step cache for EFTCAMB.
!! To perform EFTCAMB calculations, at any timestep, we need several quantities.
!! Densities, values of the EFT functions and so on.
!! Instead of complicating the interface of the functions we shall define a type that
!! holds all these informations.
!! We also need a cache for cosmological parameters that are used in the calculations.


!----------------------------------------------------------------------------------------
!> This module contains the definition of a time step cache for EFTCAMB.
!! To perform EFTCAMB calculations, at any timestep, we need several quantities.
!! Densities, values of the EFT functions and so on.
!! Instead of complicating the interface of the functions we shall define a type that
!! holds all these informations.
!! We also need a cache for cosmological parameters that are used in the calculations.

!> @author Bin Hu, Marco Raveri

module EFTCAMB_cache

    use precision
    use IniFile

    implicit none

    !----------------------------------------------------------------------------------------
    !> This is the type that defines the EFTCAMB parameter cache. The idea is to copy in here
    !! the cosmological parameters that we need from CAMB and then use this for the interfaces.
    type :: EFTCAMB_parameter_cache

        real(dl) :: omegac         !< the value of \f$ \Omega_{\rm CDM}^0 \f$.
        real(dl) :: omegab         !< the value of \f$ \Omega_{\rm b}^0 \f$.
        real(dl) :: omegav         !< the value of \f$ \Omega_{\Lambda}^0 \f$.
        real(dl) :: omegak         !< the value of \f$ \Omega_{\rm K}^0 \f$.

        real(dl) :: h0             !< reduced Hubble constant \f$ H_0/100 \f$
        real(dl) :: h0_Mpc         !< the Hubble constant in MegaParsec \f$ 10^3 \cdot H_0/c \f$

        real(dl) :: grhog          !< the value of \f$ 8 \pi G_{N} \rho_{\gamma}(t_0) \f$.
        real(dl) :: grhornomass    !< the value of \f$ 8 \pi G_{N} \rho_{\nu}(t_0) \f$.
        real(dl) :: grhormass      !< the value of \f$ 8 \pi G_{N} \rho_{m\nu}(t_0) \f$.
        real(dl) :: grhoc          !< the value of \f$ 8 \pi G_{N} \rho_{\rm CDM}(t_0) \f$.
        real(dl) :: grhob          !< the value of \f$ 8 \pi G_{N} \rho_{\rm b}(t_0) \f$.
        real(dl) :: grhov          !< the value of \f$ 8 \pi G_{N} \rho_{\Lambda}(t_0) \f$.
        real(dl) :: grhok          !< the value of \f$ 8 \pi G_{N} \rho_{\rm K}(t_0) \f$.

    contains

        procedure :: initialize => EFTCAMBParameterCacheInit  !< subroutine that initializes to zero all the elements of the parameter cache.
        procedure :: print      => EFTCAMBParameterCachePrint !< subroutine that prints the EFTCAMB parameters cache to screen.

    end type EFTCAMB_parameter_cache

    !----------------------------------------------------------------------------------------
    !> This is the type that defines the EFTCAMB time step cache.
    type :: EFTCAMB_timestep_cache

        ! 1) time:
        real(dl) :: a             !< the value of the scale factor at which the cache is being used.
        real(dl) :: tau           !< the value of conformal time at the given scale factor.
        ! 2) expansion history:
        real(dl) :: adotoa        !< the value of \f$ \mathcal{H} \f$ at the given scale factor.
        real(dl) :: Hdot          !< the value of \f$ d\mathcal{H} /d \tau \f$ at the given scale factor.
        real(dl) :: Hdotdot       !< the value of \f$ d^2 \mathcal{H} / d \tau^2 \f$ at the given scale factor.
        ! 3) matter densities:
        real(dl) :: grhoa2        !< the input value of \f$ \sum_m\rho_m / a^2 m_0^2 \f$.
        real(dl) :: grhom_t       !< the value of \f$ \sum_m\rho_m a^2 /m_0^2 \f$.
        real(dl) :: gpresm_t      !< the value of \f$ \sum_m P_m a^2 /m_0^2 \f$.
        real(dl) :: gpresdotm_t   !< the value of \f$ \sum_m\dot{P}_m a^2 /m_0^2 \f$.

        ! 4) EFT functions:
        real(dl) :: EFTOmegaV     !< the value of Omega \f$ \Omega(a) \f$.
        real(dl) :: EFTOmegaP     !< the value of the derivative wrt scale factor of Omega \f$ d \Omega(a) / da \f$.
        real(dl) :: EFTOmegaPP    !< the value of the second derivative wrt scale factor of Omega \f$ d^2 \Omega(a) / da^2 \f$.
        real(dl) :: EFTOmegaPPP   !< the value of the third derivative wrt scale factor of Omega \f$ d^3 \Omega(a) / da^3 \f$.
        real(dl) :: EFTc          !< the value of \f$ c a^2/m_0^2 \f$.
        real(dl) :: EFTcdot       !< the value of \f$ \dot{c} a^2/m_0^2 \f$.
        real(dl) :: EFTLambda     !< the value of \f$ \Lambda a^2/m_0^2 \f$.
        real(dl) :: EFTLambdadot  !< the value of \f$ \dot{\Lambda}a^2/m_0^2 \f$. Derivative of \f$ \Lambda\f$ wrt conformal time.
        real(dl) :: EFTGamma1V    !< the value of Gamma 1 \f$ \gamma_1(a) \f$.
        real(dl) :: EFTGamma1P    !< the value of the derivative wrt scale factor of Gamma 1 \f$  d \gamma_1(a) / da \f$.
        real(dl) :: EFTGamma2V    !< the value of Gamma 2 \f$ \gamma_2(a) \f$.
        real(dl) :: EFTGamma2P    !< the value of the derivative wrt scale factor of Gamma 2 \f$  d \gamma_2(a) / da \f$.
        real(dl) :: EFTGamma3V    !< the value of Gamma 3 \f$ \gamma_3(a) \f$.
        real(dl) :: EFTGamma3P    !< the value of the derivative wrt scale factor of Gamma 3 \f$  d \gamma_3(a) / da \f$.
        real(dl) :: EFTGamma4V    !< the value of Gamma 4 \f$ \gamma_4(a) \f$.
        real(dl) :: EFTGamma4P    !< the value of the derivative wrt scale factor of Gamma 4 \f$  d \gamma_4(a) / da \f$.
        real(dl) :: EFTGamma4PP   !< the value of the second derivative wrt scale factor of Gamma 4 \f$  d^2 \gamma_4(a) / da^2 \f$.
        real(dl) :: EFTGamma5V    !< the value of Gamma 5 \f$ \gamma_5(a) \f$.
        real(dl) :: EFTGamma5P    !< the value of the derivative wrt scale factor of Gamma 5 \f$  d \gamma_5(a) / da \f$.
        real(dl) :: EFTGamma6V    !< the value of Gamma 6 \f$ \gamma_6(a) \f$.
        real(dl) :: EFTGamma6P    !< the value of the derivative wrt scale factor of Gamma 6 \f$  d \gamma_6(a) / da \f$.

    contains

        procedure :: initialize => EFTCAMBTimestepCacheInit  !< subroutine that initializes to zero all the elements of the cache.

    end type EFTCAMB_timestep_cache

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes to zero all the elements of the cache.
    subroutine EFTCAMBTimestepCacheInit( self )

        implicit none

        class(EFTCAMB_timestep_cache)  :: self !< the base class.

        ! initialize all class members to zero:
        self%a            = 0._dl
        self%tau          = 0._dl
        self%adotoa       = 0._dl
        self%Hdot         = 0._dl
        self%Hdotdot      = 0._dl
        self%grhom_t      = 0._dl
        self%gpresm_t     = 0._dl
        self%gpresdotm_t  = 0._dl
        self%EFTOmegaV    = 0._dl
        self%EFTOmegaP    = 0._dl
        self%EFTOmegaPP   = 0._dl
        self%EFTOmegaPPP  = 0._dl
        self%EFTc         = 0._dl
        self%EFTcdot      = 0._dl
        self%EFTLambda    = 0._dl
        self%EFTLambdadot = 0._dl
        self%EFTGamma1V   = 0._dl
        self%EFTGamma1P   = 0._dl
        self%EFTGamma2V   = 0._dl
        self%EFTGamma2P   = 0._dl
        self%EFTGamma3V   = 0._dl
        self%EFTGamma3P   = 0._dl
        self%EFTGamma4V   = 0._dl
        self%EFTGamma4P   = 0._dl
        self%EFTGamma4PP  = 0._dl
        self%EFTGamma5V   = 0._dl
        self%EFTGamma5P   = 0._dl
        self%EFTGamma6V   = 0._dl
        self%EFTGamma6P   = 0._dl

    end subroutine EFTCAMBTimestepCacheInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes to zero all the elements of the parameter cache.
    subroutine EFTCAMBParameterCacheInit( self )

        implicit none

        class(EFTCAMB_parameter_cache)  :: self !< the base class.

        ! initialize all class members to zero:
        self%omegac      = 0._dl
        self%omegab      = 0._dl
        self%omegav      = 0._dl
        self%omegak      = 0._dl
        self%h0          = 0._dl
        self%h0_Mpc      = 0._dl
        self%grhog       = 0._dl
        self%grhornomass = 0._dl
        self%grhormass   = 0._dl
        self%grhoc       = 0._dl
        self%grhob       = 0._dl
        self%grhov       = 0._dl
        self%grhok       = 0._dl

    end subroutine EFTCAMBParameterCacheInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints the EFTCAMB parameters cache to screen.
    subroutine EFTCAMBParameterCachePrint( self )

        implicit none

        class(EFTCAMB_parameter_cache)  :: self !< the base class.

        ! print to screen the parameter cache:
        write(*,*) 'EFTCAMB parameters cache content:'
        write(*,*) 'Omega_CDM  :', self%omegac
        write(*,*) 'Omega_b    :', self%omegab
        write(*,*) 'Omega_v    :', self%omegav
        write(*,*) 'Omega_k    :', self%omegak
        write(*,*) 'h          :', self%h0
        write(*,*) 'h_Mpc      :', self%h0_Mpc
        write(*,*) 'grhog      :', self%grhog
        write(*,*) 'grnonomass :', self%grhornomass
        write(*,*) 'grhormass  :', self%grhormass
        write(*,*) 'grhoc      :', self%grhoc
        write(*,*) 'grhob      :', self%grhob
        write(*,*) 'grhov      :', self%grhov
        write(*,*) 'grhok      :', self%grhok

    end subroutine EFTCAMBParameterCachePrint

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_cache

!----------------------------------------------------------------------------------------

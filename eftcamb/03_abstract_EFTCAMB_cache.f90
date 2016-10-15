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

!> @file 03_abstract_EFTCAMB_cache.f90
!! This file contains the definition of the EFTCAMB caches.
!! These are used to store parameters that can be used by EFTCAMB, in EFTCAMB_parameter_cache
!! or are used to store partial results when solving the time evolution of perturbations,
!! in EFTCAMB_timestep_cache.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the EFTCAMB caches.
!! These are used to store parameters that can be used by EFTCAMB, in EFTCAMB_parameter_cache
!! or are used to store partial results when solving the time evolution of perturbations,
!! in EFTCAMB_timestep_cache.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_cache

    use precision
    use IniFile

    implicit none

    private

    public EFTCAMB_parameter_cache, EFTCAMB_timestep_cache

    !----------------------------------------------------------------------------------------
    !> This is the type that defines the EFTCAMB parameter cache. The idea is to copy in here
    !! the cosmological parameters that we need from CAMB and then use this for the interfaces.
    !! This also contains a wrapper for Massive Neutrinos quantities that can be used in
    !! the background code.
    type :: EFTCAMB_parameter_cache

        ! 1) relative densities:
        real(dl) :: omegac         !< the value of \f$ \Omega_{\rm CDM}^0 \f$.
        real(dl) :: omegab         !< the value of \f$ \Omega_{\rm b}^0 \f$.
        real(dl) :: omegav         !< the value of \f$ \Omega_{\Lambda}^0 \f$.
        real(dl) :: omegak         !< the value of \f$ \Omega_{\rm K}^0 \f$.
        real(dl) :: omegan         !< the value of \f$ \Omega_{m\nu}^0 \f$ of massive neutrinos.
        real(dl) :: omegag         !< the value of \f$ \Omega_{ \gamma }^0 \f$.
        real(dl) :: omegar         !< the value of \f$ \Omega_{ \nu }^0 \f$ of massless neutrinos.
        ! 2) Hubble constant:
        real(dl) :: h0             !< reduced Hubble constant \f$ H_0/100 \f$.
        real(dl) :: h0_Mpc         !< the Hubble constant in MegaParsec \f$ 10^3 \cdot H_0/c \f$.
        ! 3) densities:
        real(dl) :: grhog          !< the value of \f$ 8 \pi G_{N} \rho_{\gamma}(t_0) \f$.
        real(dl) :: grhornomass    !< the value of \f$ 8 \pi G_{N} \rho_{\nu}(t_0) \f$.
        real(dl) :: grhoc          !< the value of \f$ 8 \pi G_{N} \rho_{\rm CDM}(t_0) \f$.
        real(dl) :: grhob          !< the value of \f$ 8 \pi G_{N} \rho_{\rm b}(t_0) \f$.
        real(dl) :: grhov          !< the value of \f$ 8 \pi G_{N} \rho_{\Lambda}(t_0) \f$.
        real(dl) :: grhok          !< the value of \f$ 8 \pi G_{N} \rho_{\rm K}(t_0) \f$.
        ! 4) massive neutrinos:
        integer  :: Num_Nu_Massive                        !< number of massive neutrinos
        integer  :: Nu_mass_eigenstates                   !< number of mass eigenstates
        real(dl), allocatable, dimension(:) :: grhormass  !< densities of neutrinos in each mass eigenstate \f$ 8 \pi G_{N} \rho_{ m\nu }(t_0) \f$
        real(dl), allocatable, dimension(:) :: nu_masses  !< neutrino masses
        ! 5) massive neutrinos wrapper:
        procedure( Nu_background_Wrapper ), pointer, nopass :: Nu_background => null()  !< wrapper to the subroutine that computes the background massive neutrinos density and pressure.
        procedure( Nu_rho_Wrapper        ), pointer, nopass :: Nu_rho        => null()  !< wrapper to the subroutine that computes the background massive neutrinos density.
        procedure( Nu_pidot_Wrapper      ), pointer, nopass :: Nu_pidot      => null()  !< wrapper to the function that computes the background massive neutrinos time derivative of pressure.
        procedure( Nu_pidotdot_Wrapper   ), pointer, nopass :: Nu_pidotdot   => null()  !< wrapper to the function that computes the background massive neutrinos second time derivative of pressure.

    contains

        procedure :: initialize => EFTCAMBParameterCacheInit  !< subroutine that initializes to zero all the elements of the parameter cache.
        procedure :: print      => EFTCAMBParameterCachePrint !< subroutine that prints the EFTCAMB parameters cache to screen.
        procedure :: is_nan     => EFTCAMBParameterCacheIsNan !< Subroutine that check if an element of the EFTCAMB_parameter_cache is Nan.

    end type EFTCAMB_parameter_cache

    !----------------------------------------------------------------------------------------
    ! Interface containing the wrapper to massive neutrinos stuff.
    interface
        !----------------------------------------------------------------------------------------
        !> Wrapper to the subroutine that computes the background massive neutrinos density
        !! and pressure.
        subroutine Nu_background_Wrapper( am, rhonu, pnu )
            use precision
            implicit none
            real(dl), intent(in)  :: am     !< input scale factor times the neutrino mass
            real(dl), intent(out) :: rhonu  !< output neutrino density \f$ \frac{\rho_{\nu} a^2}{m_0^2} \f$
            real(dl), intent(out) :: pnu    !< output neutrino pressure \f$ \frac{ P_{\nu} a^2}{m_0^2} \f$
        end subroutine Nu_background_Wrapper
        !----------------------------------------------------------------------------------------
        !> Wrapper to the subroutine that computes the background massive neutrinos density.
        subroutine Nu_rho_Wrapper( am, rhonu )
            use precision
            implicit none
            real(dl), intent(in)  :: am     !< input scale factor times the neutrino mass
            real(dl), intent(out) :: rhonu  !< output neutrino density \f$ \frac{\rho_{\nu} a^2}{m_0^2} \f$
        end subroutine Nu_rho_Wrapper
        !----------------------------------------------------------------------------------------
        !> Wrapper to the function that computes the background massive neutrinos
        !! time derivative of pressure.
        function Nu_pidot_Wrapper( am, adotoa, presnu )
            use precision
            implicit none
            real(dl), intent(in)  :: am     !< input scale factor times the neutrino mass
            real(dl)              :: adotoa !< input conformal Hubble
            real(dl)              :: presnu !< input neutrino pressure
            real(dl)              :: Nu_pidot_Wrapper !< output value of the time derivative of neutrino pressure
        end function Nu_pidot_Wrapper
        !----------------------------------------------------------------------------------------
        !> Wwrapper to the function that computes the background massive neutrinos
        !! second time derivative of pressure.
        function Nu_pidotdot_Wrapper( am, adotoa, Hdot, presnu, presnudot )
            use precision
            implicit none
            real(dl), intent(in)  :: am        !< input scale factor times the neutrino mass
            real(dl)              :: adotoa    !< input conformal Hubble
            real(dl)              :: Hdot      !< input time derivative of conformal Hubble
            real(dl)              :: presnu    !< input neutrino pressure
            real(dl)              :: presnudot !< input time derivative of neutrino pressure
            real(dl)              :: Nu_pidotdot_Wrapper !< output value of the second time derivative of neutrino pressure
        end function Nu_pidotdot_Wrapper
        !----------------------------------------------------------------------------------------
    end interface


    !----------------------------------------------------------------------------------------
    !> This is the type that defines the EFTCAMB time step cache.
    type :: EFTCAMB_timestep_cache

        ! 1) time and k:
        real(dl) :: a             !< the value of the scale factor at which the cache is being used.
        real(dl) :: tau           !< the value of conformal time at the given scale factor.
        real(dl) :: k             !< the scale that is being solved for. In \f$ \mbox{Mpc}^{-1} \f$.
        ! 2) total matter densities:
        real(dl) :: grhoa2        !< the input value of \f$ \sum_m\rho_m / a^2 m_0^2 \f$.
        real(dl) :: grhom_t       !< the value of \f$ \sum_m\rho_m a^2 /m_0^2 \f$.
        real(dl) :: gpresm_t      !< the value of \f$ \sum_m P_m a^2 /m_0^2 \f$.
        real(dl) :: gpresdotm_t   !< the value of \f$ \sum_m\dot{P}_m a^2 /m_0^2 \f$.
        ! 3) densities and pressure of the various species:
        real(dl) :: grhob_t       !< the value of \f$ \rho_b a^2 / m_0^2 \f$
        real(dl) :: grhoc_t       !< the value of \f$ \rho_{cdm} a^2 / m_0^2 \f$
        real(dl) :: grhor_t       !< the value of \f$ \rho_{\nu} a^2 / m_0^2 \f$
        real(dl) :: grhog_t       !< the value of \f$ \rho_{\gamma} a^2 / m_0^2 \f$
        real(dl) :: grhov_t       !< the value of \f$ \rho_{\Lambda} a^2 / m_0^2 \f$. Used if needed, especially by designer models.
        real(dl) :: gpiv_t        !< the value of \f$ \P_{\Lambda} a^2 / m_0^2 \f$. Used if needed, especially by designer models.
        real(dl) :: grhonu_tot    !< the value of \f$ \sum_\nu \rho_{m\nu} a^2 / m_0^2 \f$
        real(dl) :: gpinu_tot     !< the value of \f$ \sum_\nu P_{m\nu} a^2 / m_0^2 \f$
        real(dl) :: grhonudot_tot !< the value of \f$ \sum_\nu \dot{\rho}_{m\nu} a^2 / m_0^2 \f$
        real(dl) :: gpinudot_tot  !< the value of \f$ \sum_\nu \dot{P}_{m\nu} a^2 / m_0^2 \f$
        ! 4) expansion history:
        real(dl) :: adotoa        !< the value of \f$ \mathcal{H} \f$ at the given scale factor.
        real(dl) :: Hdot          !< the value of \f$ d\mathcal{H} /d \tau \f$ at the given scale factor.
        real(dl) :: Hdotdot       !< the value of \f$ d^2 \mathcal{H} / d \tau^2 \f$ at the given scale factor.
        ! 5) EFT functions:
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
        ! 6) other background quantities:
        real(dl) :: grhoq         !< the value of the effective density of the Q field. Refer to the Numerical Notes for the definition.
        real(dl) :: gpresq        !< the value of the effective pressure of the Q field. Refer to the Numerical Notes for the definition.
        real(dl) :: grhodotq      !< the value of the time derivative of the effective density of the Q field. Refer to the Numerical Notes for the definition.
        real(dl) :: gpresdotq     !< the value of the time derivative of the effective pressure of the Q field. Refer to the Numerical Notes for the definition.
        ! 7) the Einstein equations coefficients:
        real(dl) :: EFTeomF       !< the value of the Einstein equations coefficient F. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomN       !< the value of the Einstein equations coefficient N. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomNdot    !< the value of the Einstein equations coefficient dN/dtau. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomX       !< the value of the Einstein equations coefficient X. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomXdot    !< the value of the Einstein equations coefficient dX/dtau. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomY       !< the value of the Einstein equations coefficient Y. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomG       !< the value of the Einstein equations coefficient G. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomU       !< the value of the Einstein equations coefficient U. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomL       !< the value of the Einstein equations coefficient L. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomM       !< the value of the Einstein equations coefficient M. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomV       !< the value of the Einstein equations coefficient V. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomVdot    !< the value of the Einstein equations coefficient dV/dtau. Refer to the Numerical Notes for the definition.
        ! 8) pi field factors:
        real(dl) :: EFTpiA1       !< the value of the pi field equation coefficient A1. Scale independent part of A. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiA2       !< the value of the pi field equation coefficient A2. Part proportional to \f$ k^2 \f$. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiB1       !< the value of the pi field equation coefficient B1. Scale independent part of B. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiB2       !< the value of the pi field equation coefficient B2. Part proportional to \f$ k^2 \f$. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiC        !< the value of the pi field equation coefficient C. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiD1       !< the value of the pi field equation coefficient D1. Scale independent part of D. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiD2       !< the value of the pi field equation coefficient D2. Part proportional to \f$ k^2 \f$. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiE        !< the value of the pi field equation coefficient E. Refer to the Numerical Notes for the definition.
        ! 9) pi field quantities:
        real(dl) :: pi            !< the value of the pi field at a given time and scale.
        real(dl) :: pidot         !< the value of the (conformal) time derivative of the pi field at a given time and scale.
        real(dl) :: pidotdot      !< the value of the (conformal) second time derivative of the pi field at a given time and scale.
        ! 10) scalar perturbations quantities:
        real(dl) :: z             !< Syncronous gauge Z perturbation.
        real(dl) :: dz            !< Syncronous gauge dot Z perturbation. This is used to store the non-RSA result.
        real(dl) :: sigma         !< Syncronous gauge sigma perturbation. This is used to store the non-RSA result.
        real(dl) :: sigmadot      !< Syncronous gauge dot sigma perturbation. This is used to store the non-RSA result.
        real(dl) :: clxg          !< Syncronous gauge radiation density perturbation.
        real(dl) :: clxr          !< Syncronous gauge massless neutrinos density perturbation.
        real(dl) :: dgpnu         !< Syncronous gauge massive neutrinos pressure perturbation.
        real(dl) :: dgrho         !< Syncronous gauge total density perturbation.
        real(dl) :: dgq           !< Syncronous gauge total velocity perturbation.
        ! 11) tensor perturbations quantities:
        real(dl) :: EFTAT         !< the value of the tensor equation coefficient A. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTBT         !< the value of the tensor equation coefficient B. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTDT         !< the value of the tensor equation coefficient D. Refer to the Numerical Notes for the definition.
        ! 12) Kinetic and Gradient quantities for the stability check:
        real(dl) :: EFT_kinetic   !< the value of the kinetic term. Refer to the Numerical Notes for the definition.
        real(dl) :: EFT_gradient  !< the value of the gradient term. Refer to the Numerical Notes for the definition.

    contains

        procedure :: initialize => EFTCAMBTimestepCacheInit  !< subroutine that initializes to zero all the elements of the cache.
        procedure :: is_nan     => EFTCAMBTimestepCacheIsNan !< Subroutine that check if an element of the EFTCAMB_timestep_cache is Nan.

    end type EFTCAMB_timestep_cache

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes to zero all the elements of the cache.
    subroutine EFTCAMBTimestepCacheInit( self )

        implicit none

        class(EFTCAMB_timestep_cache)  :: self !< the base class.

        ! initialize all class members to zero:
        ! 1) time and k:
        self%a             = 0._dl
        self%tau           = 0._dl
        self%k             = 0._dl
        ! 2) expansion history:
        self%adotoa        = 0._dl
        self%Hdot          = 0._dl
        self%Hdotdot       = 0._dl
        ! 3) total matter densities:
        self%grhoa2        = 0._dl
        self%grhom_t       = 0._dl
        self%gpresm_t      = 0._dl
        self%gpresdotm_t   = 0._dl
        ! 4) densities and pressure of the various species:
        self%grhob_t       = 0._dl
        self%grhoc_t       = 0._dl
        self%grhor_t       = 0._dl
        self%grhog_t       = 0._dl
        self%grhov_t       = 0._dl
        self%gpiv_t        = 0._dl
        self%grhonu_tot    = 0._dl
        self%gpinu_tot     = 0._dl
        self%grhonudot_tot = 0._dl
        self%gpinudot_tot  = 0._dl
        ! 5) EFT functions:
        self%EFTOmegaV     = 0._dl
        self%EFTOmegaP     = 0._dl
        self%EFTOmegaPP    = 0._dl
        self%EFTOmegaPPP   = 0._dl
        self%EFTc          = 0._dl
        self%EFTcdot       = 0._dl
        self%EFTLambda     = 0._dl
        self%EFTLambdadot  = 0._dl
        self%EFTGamma1V    = 0._dl
        self%EFTGamma1P    = 0._dl
        self%EFTGamma2V    = 0._dl
        self%EFTGamma2P    = 0._dl
        self%EFTGamma3V    = 0._dl
        self%EFTGamma3P    = 0._dl
        self%EFTGamma4V    = 0._dl
        self%EFTGamma4P    = 0._dl
        self%EFTGamma4PP   = 0._dl
        self%EFTGamma5V    = 0._dl
        self%EFTGamma5P    = 0._dl
        self%EFTGamma6V    = 0._dl
        self%EFTGamma6P    = 0._dl
        ! 6) other background quantities:
        self%grhoq         = 0._dl
        self%gpresq        = 0._dl
        self%grhodotq      = 0._dl
        self%gpresdotq     = 0._dl
        ! 7) the Einstein equations coefficients:
        self%EFTeomF       = 0._dl
        self%EFTeomN       = 0._dl
        self%EFTeomNdot    = 0._dl
        self%EFTeomX       = 0._dl
        self%EFTeomXdot    = 0._dl
        self%EFTeomY       = 0._dl
        self%EFTeomG       = 0._dl
        self%EFTeomU       = 0._dl
        self%EFTeomL       = 0._dl
        self%EFTeomM       = 0._dl
        self%EFTeomV       = 0._dl
        self%EFTeomVdot    = 0._dl
        ! 8) pi field factors:
        self%EFTpiA1       = 0._dl
        self%EFTpiA2       = 0._dl
        self%EFTpiB1       = 0._dl
        self%EFTpiB2       = 0._dl
        self%EFTpiC        = 0._dl
        self%EFTpiD1       = 0._dl
        self%EFTpiD2       = 0._dl
        self%EFTpiE        = 0._dl
        ! 9) pi field quantities:
        self%pi            = 0._dl
        self%pidot         = 0._dl
        self%pidotdot      = 0._dl
        ! 10) perturbations quantities:
        self%z             = 0._dl
        self%clxg          = 0._dl
        self%clxr          = 0._dl
        self%dgpnu         = 0._dl
        self%dgrho         = 0._dl
        self%dgq           = 0._dl
        ! 11) tensor perturbations quantities:
        self%EFTAT         = 0._dl
        self%EFTBT         = 0._dl
        self%EFTDT         = 0._dl
        ! 12) Kinetic and Gradient quantities for the stability check:
        self%EFT_kinetic   = 0._dl
        self%EFT_gradient  = 0._dl

    end subroutine EFTCAMBTimestepCacheInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that check if an element of the EFTCAMB_timestep_cache is Nan.
    subroutine EFTCAMBTimestepCacheIsNan( self, HaveNan )

      implicit none

      class(EFTCAMB_timestep_cache), intent(in)  :: self    !< The base class.
      logical, intent(out)                       :: HaveNan !< Logical variable which describes the presence of a Nan variable.
                                                            !< If an element of the EFTCAMB_timestep_cache is Nan, you get HaveNan=.True.

      HaveNan = .False.
      HaveNan = HaveNan.or.IsNaN(self%a)
      HaveNan = HaveNan.or.IsNaN(self%tau)
      HaveNan = HaveNan.or.IsNaN(self%k)
      HaveNan = HaveNan.or.IsNaN(self%adotoa)
      HaveNan = HaveNan.or.IsNaN(self%Hdot)
      HaveNan = HaveNan.or.IsNaN(self%Hdotdot)
      HaveNan = HaveNan.or.IsNaN(self%grhoa2)
      HaveNan = HaveNan.or.IsNaN(self%grhom_t)
      HaveNan = HaveNan.or.IsNaN(self%gpresm_t)
      HaveNan = HaveNan.or.IsNaN(self%gpresdotm_t)
      HaveNan = HaveNan.or.IsNaN(self%grhob_t)
      HaveNan = HaveNan.or.IsNaN(self%grhoc_t)
      HaveNan = HaveNan.or.IsNaN(self%grhor_t)
      HaveNan = HaveNan.or.IsNaN(self%grhog_t)
      HaveNan = HaveNan.or.IsNaN(self%grhov_t)
      HaveNan = HaveNan.or.IsNaN(self%gpiv_t)
      HaveNan = HaveNan.or.IsNaN(self%grhonu_tot)
      HaveNan = HaveNan.or.IsNaN(self%gpinu_tot)
      HaveNan = HaveNan.or.IsNaN(self%grhonudot_tot)
      HaveNan = HaveNan.or.IsNaN(self%gpinudot_tot)
      HaveNan = HaveNan.or.IsNaN(self%EFTOmegaV)
      HaveNan = HaveNan.or.IsNaN(self%EFTOmegaP)
      HaveNan = HaveNan.or.IsNaN(self%EFTOmegaPP)
      HaveNan = HaveNan.or.IsNaN(self%EFTOmegaPPP)
      HaveNan = HaveNan.or.IsNaN(self%EFTc)
      HaveNan = HaveNan.or.IsNaN(self%EFTcdot)
      HaveNan = HaveNan.or.IsNaN(self%EFTLambda)
      HaveNan = HaveNan.or.IsNaN(self%EFTLambdadot)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma1V)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma1P)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma2V)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma2P)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma3V)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma3P)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma4V)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma4P)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma4PP)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma5V)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma5P)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma6V)
      HaveNan = HaveNan.or.IsNaN(self%EFTGamma6P)
      HaveNan = HaveNan.or.IsNaN(self%grhoq)
      HaveNan = HaveNan.or.IsNaN(self%gpresq)
      HaveNan = HaveNan.or.IsNaN(self%grhodotq)
      HaveNan = HaveNan.or.IsNaN(self%gpresdotq)
      HaveNan = HaveNan.or.IsNaN(self%EFTeomF)
      HaveNan = HaveNan.or.IsNaN(self%EFTeomN)
      HaveNan = HaveNan.or.IsNaN(self%EFTeomNdot)
      HaveNan = HaveNan.or.IsNaN(self%EFTeomX)
      HaveNan = HaveNan.or.IsNaN(self%EFTeomXdot)
      HaveNan = HaveNan.or.IsNaN(self%EFTeomY)
      HaveNan = HaveNan.or.IsNaN(self%EFTeomG)
      HaveNan = HaveNan.or.IsNaN(self%EFTeomU)
      HaveNan = HaveNan.or.IsNaN(self%EFTeomL)
      HaveNan = HaveNan.or.IsNaN(self%EFTeomM)
      HaveNan = HaveNan.or.IsNaN(self%EFTeomV)
      HaveNan = HaveNan.or.IsNaN(self%EFTeomVdot)
      HaveNan = HaveNan.or.IsNaN(self%EFTpiA1)
      HaveNan = HaveNan.or.IsNaN(self%EFTpiA2)
      HaveNan = HaveNan.or.IsNaN(self%EFTpiB1)
      HaveNan = HaveNan.or.IsNaN(self%EFTpiB2)
      HaveNan = HaveNan.or.IsNaN(self%EFTpiC)
      HaveNan = HaveNan.or.IsNaN(self%EFTpiD1)
      HaveNan = HaveNan.or.IsNaN(self%EFTpiD2)
      HaveNan = HaveNan.or.IsNaN(self%EFTpiE)
      HaveNan = HaveNan.or.IsNaN(self%pi)
      HaveNan = HaveNan.or.IsNaN(self%pidot)
      HaveNan = HaveNan.or.IsNaN(self%pidotdot)
      HaveNan = HaveNan.or.IsNaN(self%z)
      HaveNan = HaveNan.or.IsNaN(self%clxg)
      HaveNan = HaveNan.or.IsNaN(self%clxr)
      HaveNan = HaveNan.or.IsNaN(self%dgpnu)
      HaveNan = HaveNan.or.IsNaN(self%dgrho)
      HaveNan = HaveNan.or.IsNaN(self%dgq)
      HaveNan = HaveNan.or.IsNaN(self%EFTAT)
      HaveNan = HaveNan.or.IsNaN(self%EFTBT)
      HaveNan = HaveNan.or.IsNaN(self%EFTDT)
      HaveNan = HaveNan.or.IsNaN(self%EFT_kinetic)
      HaveNan = HaveNan.or.IsNaN(self%EFT_gradient)

    end subroutine EFTCAMBTimestepCacheIsNan

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
        self%omegan      = 0._dl
        self%omegag      = 0._dl
        self%omegar      = 0._dl
        self%h0          = 0._dl
        self%h0_Mpc      = 0._dl
        self%grhog       = 0._dl
        self%grhornomass = 0._dl
        self%grhoc       = 0._dl
        self%grhob       = 0._dl
        self%grhov       = 0._dl
        self%grhok       = 0._dl
        self%Num_Nu_Massive       = 0
        self%Nu_mass_eigenstates  = 0
        if ( allocated(self%grhormass) ) deallocate(self%grhormass)
        if ( allocated(self%nu_masses) ) deallocate(self%nu_masses)
        if ( associated(self%Nu_background) ) nullify(self%Nu_background)
        if ( associated(self%Nu_rho)        ) nullify(self%Nu_rho)
        if ( associated(self%Nu_pidot)      ) nullify(self%Nu_pidot)
        if ( associated(self%Nu_pidotdot)   ) nullify(self%Nu_pidotdot)

    end subroutine EFTCAMBParameterCacheInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints the EFTCAMB parameters cache to screen.
    subroutine EFTCAMBParameterCachePrint( self )

        implicit none

        class(EFTCAMB_parameter_cache)  :: self !< the base class.

        integer  :: i

        ! print to screen the parameter cache:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') 'EFTCAMB parameters cache content:'
        write(*,'(a)') "***************************************************************"
        write(*,'(a14,E13.6)') ' Omega_CDM  : ', self%omegac
        write(*,'(a14,E13.6)') ' Omega_b    : ', self%omegab
        write(*,'(a14,E13.6)') ' Omega_v    : ', self%omegav
        write(*,'(a14,E13.6)') ' Omega_k    : ', self%omegak
        write(*,'(a14,E13.6)') ' Omega_n    : ', self%omegan
        write(*,'(a14,E13.6)') ' Omega_g    : ', self%omegag
        write(*,'(a14,E13.6)') ' Omega_r    : ', self%omegar
        write(*,'(a14,F12.6)') ' h          : ', self%h0
        write(*,'(a14,E13.6)') ' h_Mpc      : ', self%h0_Mpc
        write(*,'(a14,E13.6)') ' grhog      : ', self%grhog
        write(*,'(a14,E13.6)') ' grnonomass : ', self%grhornomass
        write(*,'(a14,E13.6)') ' grhoc      : ', self%grhoc
        write(*,'(a14,E13.6)') ' grhob      : ', self%grhob
        write(*,'(a14,E13.6)') ' grhov      : ', self%grhov
        write(*,'(a14,E13.6)') ' grhok      : ', self%grhok
        write(*,'(a22,I10)') ' Num_Nu_Massive      : ', self%Num_Nu_Massive
        write(*,'(a22,I10)') ' Nu_mass_eigenstates : ', self%Nu_mass_eigenstates
        do i=1, self%Nu_mass_eigenstates
            write(*,'(a11,I3,a9,E13.6)') ' grhormass(',i,')      : ', self%grhormass(i)
            write(*,'(a11,I3,a9,E13.6)') ' nu_masses(',i,')      : ', self%nu_masses(i)
        end do
        write(*,'(a)') "***************************************************************"

    end subroutine EFTCAMBParameterCachePrint

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that check if an element of the EFTCAMB_parameter_cache is Nan.
   subroutine EFTCAMBParameterCacheIsNan( self, HaveNan )

     implicit none

     class(EFTCAMB_parameter_cache), intent(in)  :: self    !< The base class.
     logical, intent(out)                        :: HaveNan !< Logical variable which describes the presence of a Nan variable.
                                                            !< If an element of the EFTCAMB_parameter_cache is Nan, you get HaveNan=.True.
     integer                                     :: i

     HaveNan = .False.
     HaveNan = HaveNan.or.IsNaN(self%omegac)
     HaveNan = HaveNan.or.IsNaN(self%omegab)
     HaveNan = HaveNan.or.IsNaN(self%omegav)
     HaveNan = HaveNan.or.IsNaN(self%omegak)
     HaveNan = HaveNan.or.IsNaN(self%omegan)
     HaveNan = HaveNan.or.IsNaN(self%omegag)
     HaveNan = HaveNan.or.IsNaN(self%omegar)
     HaveNan = HaveNan.or.IsNaN(self%h0)
     HaveNan = HaveNan.or.IsNaN(self%h0_Mpc)
     HaveNan = HaveNan.or.IsNaN(self%grhog)
     HaveNan = HaveNan.or.IsNaN(self%grhornomass)
     HaveNan = HaveNan.or.IsNaN(self%grhoc)
     HaveNan = HaveNan.or.IsNaN(self%grhob)
     HaveNan = HaveNan.or.IsNaN(self%grhov)
     HaveNan = HaveNan.or.IsNaN(self%grhok)
     HaveNan = HaveNan.or.IsNaN(self%Num_Nu_Massive*1.0_dl)
     HaveNan = HaveNan.or.IsNaN(self%Nu_mass_eigenstates*1.0)
     !
     do i=1, self%Nu_mass_eigenstates
       HaveNan = HaveNan.or.IsNaN(self%grhormass(i)).or.IsNaN(self%nu_masses(i))
     end do

   end subroutine EFTCAMBParameterCacheIsNan

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_cache

!----------------------------------------------------------------------------------------

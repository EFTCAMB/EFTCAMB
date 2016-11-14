! -------------------------------------------------------------------------------------------------
!
!   EFTCAMB
!
!   Developed and implemented by:
!       Bin Hu (hu@lorentz.leidenuniv.nl), Marco Raveri (mraveri@sissa.it)
!       Noemi Frusciante (fruscian@iap.fr), Alessandra Silvestri (silvestri@lorentz.leidenuniv.nl)
!
!
!   Compile time EFTCAMB options and Dark Energy equation of state module.
!
!   For more informations about the methods contained in this file see
!   the documentation: arXiv:1405.3590
!
! -------------------------------------------------------------------------------------------------

! Definition of all the EFTCAMB compile time flags and physical constants later used by the code.
module EFTDef

    use Precision
    use ModelParams

    implicit none
    ! EFT compile time flags:

    ! 1) Turn on pi field flag:
    !    Sets the scale factor at which the code starts to evolve the pi field.
    !    At times earlier than these the code evolves perturbations as in GR.
    !    This number is used as a lower bound and is refined if the teory is very close to GR
    !    by the EFTreturntoGR module.
    real(dl), parameter :: EFTturnonpiInitial = 1.d-2

    ! 2) Return to GR flag:
    !    This is the threshold at which a theory is considered to be exactly GR.
    real(dl), parameter :: EFTtoGR = 1.d-8

    ! 3) Early time stability:
    !    EFTCAMB checks the stability of the theory at all times making sure that the choosen model
    !    is viable during radiation, matter and DE eras.
    !    It is possible to enforce stability only at late times (from EFTturnonpiInitial to today)
    !    by setting this flag to false.
    !    This choice will however make the results dependent on what one calls late time,
    !    i.e. the choice of EFTturnonpiInitial, so we advice to state it clearly when reporting results.
    logical, parameter :: EarlyTimeStability = .true.

    ! 4) Priors flags:
    !    These flag will decide wether to impose a prior or not. For a detailed explanation refer to
    !    the documentation: arXiv:1405.3590

    !    1- Decides wether to use the old stability ghost and gradient conditions or the new ones.
    logical, parameter :: EFT_old_stability = .false.
    !    2- Speed of light prior (physical prior):
    !       the propagation speed of pi should be smaller than the speed of light.
    logical, parameter :: EFTlightspeedPrior = .false.
    !    3- Effective pi mass prior (physical prior):
    !       the effective mass of pi should be greater than zero.
    logical, parameter :: EFTpiMassPrior = .false.

    ! 5) EFTCAMB debug flags:
    !    These flags are used to deeply debug EFTCAMB.
    !    Don't turn them on unless you know what they are going to do.
    !    These flags have to be used in conjunction with two other changes.
    !           1- The fixq flag in cmbmain should be set different from zero.
    !           2- The path of the files created within cmbmain should be consistent
    !              (search EFTCAMB in cmbmain).
    !    To debug the code it should be launched on just one thread by
    !    changing the appropriate entry in the parameter file.
    logical, parameter :: DebugEFT       = .false.
    logical, parameter :: DebugEFTOutput = .false.
    logical, parameter :: DebugDesigner  = .false.
    logical, parameter :: DebugMapping   = .false.

    ! EFT cosmological and physical constants:

    ! 1) Physical constants:
    real(dl), parameter :: const_pi_EFT = 3.1415926535897932384626433832795_dl
    real(dl), parameter :: c_EFT = 2.99792458e8_dl
    real(dl), parameter :: G_EFT = 6.6738e-11_dl
    real(dl), parameter :: Mpc_EFT = 3.085678e22_dl
    real(dl), parameter :: sigma_boltz_EFT = 5.6704e-8_dl
    real(dl), parameter :: kappa_EFT =8._dl*const_pi_EFT*G_EFT

    ! EFTCAMB auxiliary flags: these are handeled by the program at runtime.

    ! 1) Effective turn on pi field flag:
    !    This is extablished by the code as the time at which the model considered
    !    is so close to GR that there is no seizable DE/MG effect. If such a time exists then
    !    it will override the previous value otherwise it will be equal to that value.
    real(dl) :: EFTturnonpi

    ! 2) Use with cosmomc flag:
    !    If EFTCAMB is used inside CosmoMC then the initialization subroutine has to be called
    !    inside it so we don't need to call it again inside CAMB.
    !    This flag takes care of it. Do not change it, as the
    !    CosmoMC code will take care of handeling it.
    integer :: EFTCAMBuseinCOSMOMC = 1

end module EFTDef

! -------------------------------------------------------------------------------------------------

! Definition of the model selection flags at compile time. This will compile EFTCAMB that will work
! for the specified model only. On the other hand this will allow the compiler to optimize more
! the code resulting in performances gain.
module compile_time_eft

    use Precision
    use ModelParams

    ! the behaviour of the module is decided by this flag. If you want to turn it on and define
    ! the model at compile time turn it to true and change the model selection flags according
    ! to your needs.
    logical, parameter :: compile_time_eftcamb = .false.

    integer, parameter :: CT_EFTflag = 0

    integer, parameter :: CT_PureEFTmodelOmega  = 0
    integer, parameter :: CT_PureEFTmodelGamma1 = 0
    integer, parameter :: CT_PureEFTmodelGamma2 = 0
    integer, parameter :: CT_PureEFTmodelGamma3 = 0
    integer, parameter :: CT_PureEFTmodelGamma4 = 0
    integer, parameter :: CT_PureEFTmodelGamma5 = 0
    integer, parameter :: CT_PureEFTmodelGamma6 = 0

    integer, parameter :: CT_DesignerEFTmodel     = 1
    integer, parameter :: CT_AltParEFTmodel       = 1
    integer, parameter :: CT_FullMappingEFTmodel  = 1

    integer, parameter :: CT_EFTwDE = 0
    integer, parameter :: CT_PureEFTHorndeski = 0
    integer, parameter :: CT_RPHmassPmodel = 0
    integer, parameter :: CT_RPHkineticitymodel = 0
    integer, parameter :: CT_RPHbraidingmodel = 0
    integer, parameter :: CT_RPHtensormodel = 0
    logical, parameter :: CT_HoravaSolarSystem = .false.

end module compile_time_eft

! -------------------------------------------------------------------------------------------------

! Definition of the dark energy equation of state (w_DE) used by the designer and pure EFT code.
module EFTdeEOS

    use EFTDef
    implicit none

contains

    !DEC$ ATTRIBUTES FORCEINLINE :: EFTw
    function EFTw(a,deriv)
        implicit none
        real(dl), intent(in) :: a
        integer , intent(in) :: deriv
        real(dl) :: EFTw, a0
        !Protection against numerical errors at very early times
        if (a .lt. 1.d-8) then
            a0 = 1.d-8
        else
            a0 = a
        end if
        ! Definition of the function itself
        !   The user has to define:
        !   deriv= 0: the function w(a) itself;
        !   deriv= 1: the first derivative dw(a)/da wrt the scale factor;
        !   deriv= 2: the second derivative d^2w(a)/da^2 wrt the scale factor;
        !   deriv= 3: the integral (see the numerical notes).
        select case (CP%EFTwDE)
            case (0) ! LCDM
                if (deriv==0) then
                    EFTw = -1._dl
                else if (deriv==1) then
                    EFTw = 0._dl
                else if (deriv==2) then
                    EFTw = 0._dl
                else if (deriv==3) then
                    EFTw = a0*a0
                end if
            case (1) ! wCDM
                if (deriv==0) then
                    EFTw = CP%EFTw0
                else if (deriv==1) then
                    EFTw = 0._dl
                else if (deriv==2) then
                    EFTw = 0._dl
                else if (deriv==3) then
                    EFTw = a0**(-1._dl-3._dl*CP%EFTw0)
                end if
            case (2) ! CPL
                if (deriv==0) then
                    EFTw = CP%EFTw0 + CP%EFTwa*(1._dl-a0)
                else if (deriv==1) then
                    EFTw = -CP%EFTwa
                else if (deriv==2) then
                    EFTw = 0._dl
                else if (deriv==3) then
                    EFTw = a0**(-1._dl-3._dl*(+CP%EFTw0+CP%EFTwa))*exp(3._dl*CP%EFTwa*(a0-1._dl))
                end if
            case (3) ! JPB
                if (deriv==0) then
                    EFTw = CP%EFTw0 +CP%EFTwa*(1._dl-a0)*a0**(CP%EFTwn-1._dl)
                else if (deriv==1) then
                    EFTw = -CP%EFTwa*a0**(CP%EFTwn-2._dl)*(1._dl+(a0-1._dl)*CP%EFTwn)
                else if (deriv==2) then
                    EFTw = -CP%EFTwa*a0**(CP%EFTwn-3._dl)*(CP%EFTwn-1._dl)*(2._dl+(a0-1._dl)*CP%EFTwn)
                else if (deriv==3) then
                    if (CP%EFTwn/=1) then
                        EFTw = a0**(-1._dl-3._dl*CP%EFTw0)*exp((3._dl*(a0+a0**CP%EFTwn*(a0*(CP%EFTwn-1._dl)&
                            &-CP%EFTwn))*CP%EFTwa)/(a0*(CP%EFTwn-1._dl)*CP%EFTwn))
                    else if (CP%EFTwn==1) then
                        EFTw = a0**(-1._dl -3._dl*(+CP%EFTw0+CP%EFTwa))*exp(3._dl*CP%EFTwa*(a0-1._dl))
                    end if
                end if
            case (4) ! Inflection model
                if (deriv==0) then
                    EFTw = CP%EFTw0 +CP%EFTwa*(CP%EFTwat-a0)**2
                else if (deriv==1) then
                    EFTw = 2._dl*CP%EFTwa*(a0-CP%EFTwat)
                else if (deriv==2) then
                    EFTw = 2._dl*CP%EFTwa
                else if (deriv==3) then
                    EFTw =a0**(-1._dl-3._dl*(CP%EFTw0+CP%EFTwat**2*CP%EFTwa))*&
                        &exp(-1.5_dl*(a0-1._dl)*(1._dl+a0-4._dl*CP%EFTwat)*CP%EFTwa)
                end if
            case (5) ! Taylor expansion
                if (deriv==0) then
                    EFTw = CP%EFTw0 +CP%EFTwa*a0 +0.5_dl*a0**2*CP%EFTw2 + a0**3*CP%EFTw3/6._dl
                else if (deriv==1) then
                    EFTw = CP%EFTwa +CP%EFTw2*a0 +0.5_dl*a0**2*CP%EFTw3
                else if (deriv==2) then
                    EFTw = CP%EFTw2 +CP%EFTw3*a0
                else if (deriv==3) then
                    EFTw = a0**(-1._dl -3._dl*CP%EFTw0)*exp(-(a0-1._dl)/12._dl*(9._dl*(1._dl+a0)&
                        &*CP%EFTw2+2._dl*(1._dl+a0+a0**2)*CP%EFTw3 +36._dl*CP%EFTwa))
                end if
            case (6) ! User defined
                write(*,*) 'EFTCAMB: you have to define EFTw. Go to EFTdeEOS in EFT_def.f90'
                stop
            case default
                write(*,*) 'EFTCAMB: you have to select a model for the dark energy EOS.'
                stop
        end select
    end function EFTw

end module EFTdeEOS

! -------------------------------------------------------------------------------------------------

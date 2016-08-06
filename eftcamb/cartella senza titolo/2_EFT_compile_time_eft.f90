

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

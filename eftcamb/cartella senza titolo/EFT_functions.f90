! -------------------------------------------------------------------------------------------------
!
!   EFTCAMB
!
!   Developed and implemented by:
!       Bin Hu (hu@lorentz.leidenuniv.nl), Marco Raveri (mraveri@sissa.it)
!       Noemi Frusciante (fruscian@iap.fr), Alessandra Silvestri (silvestri@lorentz.leidenuniv.nl)
!
!
!   EFTCAMB module that defines the EFT functions for the choosen model.
!
!   For more informations about the methods contained in this file see
!   the documentation: arXiv:1405.3590
!
! -------------------------------------------------------------------------------------------------


! Definitions of the EFT functions later used by the code.
! After the designer or mapping code have been called the EFT functions can be used directly from
! this module.
module EFTfunctions
    use EFTDef
    use EFTdeEOS
    use EFTdesigner
    use RPH_to_EFT

contains

    ! ----------------------------

    function EFTOmega(a,deriv)
        !EFT function Omega
        implicit none
        real(dl) :: EFTOmega
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%EFTflag)
            case (0) ! GR
                EFTOmega = 0._dl
            case (1) ! Pure EFT
                select case (CP%PureEFTmodelOmega)
                    case (0) ! 0) Operator neglected;
                        EFTOmega = 0._dl
                    case (1) ! 1) Constant model;
                        if (deriv==0) then
                            EFTOmega = CP%EFTOmega0
                        else
                            EFTOmega = 0._dl
                        end if
                    case (2) ! 2) Linear model;
                        if (deriv==0) then
                            EFTOmega= CP%EFTOmega0*a
                        else if (deriv==1) then
                            EFTOmega= CP%EFTOmega0
                        else
                            EFTOmega= 0._dl
                        end if
                    case (3) ! 3) Power law model;
                        if (deriv==0) then
                            EFTOmega= CP%EFTOmega0*a**CP%EFTOmegaExp
                        else if (deriv==1) then
                            EFTOmega= CP%EFTOmegaExp*CP%EFTOmega0*a**(CP%EFTOmegaExp-1.0_dl)
                        else if (deriv==2) then
                            EFTOmega= (CP%EFTOmegaExp-1.0_dl)*CP%EFTOmegaExp*CP%EFTOmega0*&
                                &a**(CP%EFTOmegaExp-2.0_dl)
                        else if (deriv==3) then
                            EFTOmega= (CP%EFTOmegaExp-2.0_dl)*(CP%EFTOmegaExp-1.0_dl)*&
                                &CP%EFTOmegaExp*CP%EFTOmega0*a**(CP%EFTOmegaExp-3.0_dl)
                        end if
                    case (4) ! 4) Exponential model;
                        if (deriv==0) then
                            EFTOmega= Exp(CP%EFTOmega0*a**CP%EFTOmegaExp)-1._dl
                        else if (deriv==1) then
                            EFTOmega= CP%EFTOmegaExp*CP%EFTOmega0*a**(CP%EFTOmegaExp-1.0_dl)&
                                &*Exp(CP%EFTOmega0*a**CP%EFTOmegaExp)
                        else if (deriv==2) then
                            EFTOmega= a**(-2._dl+CP%EFTOmegaExp)*Exp(a**CP%EFTOmegaExp*CP%EFTOmega0)*&
                                &CP%EFTOmega0*CP%EFTOmegaExp*(-1._dl+CP%EFTOmegaExp&
                                &+a**CP%EFTOmegaExp*CP%EFTOmega0*CP%EFTOmegaExp)
                        else if (deriv==3) then
                            EFTOmega= a**(-3._dl+CP%EFTOmegaExp)*Exp(a**CP%EFTOmegaExp*CP%EFTOmega0)*&
                                &CP%EFTOmega0*CP%EFTOmegaExp*(2._dl-3._dl*(1._dl+a**CP%EFTOmegaExp*&
                                &CP%EFTOmega0)*CP%EFTOmegaExp+(1._dl+3._dl*a**CP%EFTOmegaExp*&
                                &CP%EFTOmega0+a**(2._dl*CP%EFTOmegaExp)*CP%EFTOmega0**2)&
                                &*CP%EFTOmegaExp**2)
                        end if
                    case (5) ! 5) User defined;
                        stop 'You have to define your model for EFTOmega first! Go to the file EFT_functions.f90.'
                        if (deriv==0) then
                        else if (deriv==1) then
                        else if (deriv==2) then
                        else if (deriv==3) then
                        end if
                    case default
                        stop 'Wrong model for EFTOmega.'
                end select
            case (2) ! Designer mapping EFT
                select case (CP%DesignerEFTModel)
                    case (1) ! Designer f(R) models
                        if (deriv==0) then
                            EFTOmega= Designed_EFT_Function(OmegaDes, a, 0._dl)
                        else if (deriv==1) then
                            EFTOmega= Designed_EFT_Function(OmegapDes, a, 0._dl)
                        else if (deriv==2) then
                            EFTOmega= Designed_EFT_Function(OmegappDes, a, 0._dl)
                        else if (deriv==3) then
                            EFTOmega= Designed_EFT_Function(Omega3pDes, a, 0._dl)
                        end if
                    case (2) ! Designer mc quintessence
                        EFTOmega = 0._dl
                    case default
                        stop 'Wrong model for EFTOmega.'
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH Horndeski
                        EFTOmega = RPH_EFT_Omega(a, deriv)
                    case default
                        stop 'Wrong model for EFTOmega.'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        if (deriv==0) then
                            EFTOmega = CP%Horava_eta/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                        else
                            EFTOmega = 0._dl
                        end if
                    case default
                        stop 'Wrong model for EFTOmega.'
                end select
            case default
                stop 'Wrong model for EFTOmega.'
        end select

    end function EFTOmega

    ! ----------------------------

    function EFTGamma1(a,deriv)
        ! EFT function Gamma_1
        implicit none
        real(dl) :: EFTGamma1
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%EFTflag)
            case (0) ! GR
                EFTGamma1 = 0._dl
            case (1) ! Pure EFT
                select case (CP%PureEFTmodelGamma1)
                    case (0) ! 0) Operator neglected;
                        EFTGamma1 = 0._dl
                    case (1) ! 1) Constant model;
                        if (deriv==0) then
                            EFTGamma1 = CP%EFTGamma10
                        else
                            EFTGamma1 = 0._dl
                        end if
                    case (2) ! 2) Linear model;
                        if (deriv==0) then
                            EFTGamma1= CP%EFTGamma10*a
                        else if (deriv==1) then
                            EFTGamma1= CP%EFTGamma10
                        else
                            EFTGamma1= 0._dl
                        end if
                    case (3) ! 3) Power law model;
                        if (deriv==0) then
                            EFTGamma1= CP%EFTGamma10*a**CP%EFTGamma1Exp
                        else if (deriv==1) then
                            EFTGamma1= CP%EFTGamma1Exp*CP%EFTGamma10*a**(CP%EFTGamma1Exp-1.0_dl)
                        else if (deriv==2) then
                            EFTGamma1= (CP%EFTGamma1Exp-1.0_dl)*CP%EFTGamma1Exp*CP%EFTGamma10*&
                                &a**(CP%EFTGamma1Exp-2.0_dl)
                        else if (deriv==3) then
                            EFTGamma1= (CP%EFTGamma1Exp-2.0_dl)*(CP%EFTGamma1Exp-1.0_dl)*&
                                &CP%EFTGamma1Exp*CP%EFTGamma10*a**(CP%EFTGamma1Exp-3.0_dl)
                        end if
                    case (4) ! 4) Exponential model;
                        if (deriv==0) then
                            EFTGamma1= Exp(CP%EFTGamma10*a**CP%EFTGamma1Exp)-1._dl
                        else if (deriv==1) then
                            EFTGamma1= CP%EFTGamma1Exp*CP%EFTGamma10*a**(CP%EFTGamma1Exp-1.0_dl)&
                                &*Exp(CP%EFTGamma10*a**CP%EFTGamma1Exp)
                        else if (deriv==2) then
                            EFTGamma1= a**(-2._dl+CP%EFTGamma1Exp)*Exp(a**CP%EFTGamma1Exp*CP%EFTGamma10)*&
                                &CP%EFTGamma10*CP%EFTGamma1Exp*(-1._dl+CP%EFTGamma1Exp&
                                &+a**CP%EFTGamma1Exp*CP%EFTGamma10*CP%EFTGamma1Exp)
                        else if (deriv==3) then
                            EFTGamma1= a**(-3._dl+CP%EFTGamma1Exp)*Exp(a**CP%EFTGamma1Exp*CP%EFTGamma10)*&
                                &CP%EFTGamma10*CP%EFTGamma1Exp*(2._dl-3._dl*(1._dl+a**CP%EFTGamma1Exp*&
                                &CP%EFTGamma10)*CP%EFTGamma1Exp+(1._dl+3._dl*a**CP%EFTGamma1Exp*&
                                &CP%EFTGamma10+a**(2._dl*CP%EFTGamma1Exp)*CP%EFTGamma10**2)&
                                &*CP%EFTGamma1Exp**2)
                        end if
                    case (5) ! 5) User defined;
                        stop 'You have to define your model for EFTGamma1 first! Go to the file EFT_functions.f90.'
                        if (deriv==0) then
                        else if (deriv==1) then
                        else if (deriv==2) then
                        else if (deriv==3) then
                        end if
                    case default
                        stop 'Wrong choice for EFTGamma1.'
                end select
            case (2) ! Designer mapping EFT
                select case (CP%DesignerEFTModel)
                    case (1) ! Designer f(R) models
                        EFTGamma1=0._dl
                    case (2) ! Designer mc quintessence
                        EFTGamma1=0._dl
                    case default
                        stop 'Wrong choice for EFTGamma1.'
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH Horndeski
                        EFTGamma1 = RPH_EFT_Gamma1(a, deriv)
                    case default
                        stop 'Wrong choice for EFTGamma1.'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        EFTGamma1 = 0._dl ! for Horava EFTGamma1 is obtained with the eom in the appropriate places
                    case default
                        stop 'Wrong model for EFTGamma1.'
                end select
            case default
                stop 'Wrong choice for EFTGamma1.'
        end select

    end function EFTGamma1

    ! ----------------------------

    function EFTGamma2(a,deriv)
        ! EFT function Gamma_2
        implicit none
        real(dl) :: EFTGamma2
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%EFTflag)
            case (0) ! GR
                EFTGamma2 = 0._dl
            case (1) ! Pure EFT
                select case (CP%PureEFTmodelGamma2)
                    case (0) ! 0) Operator neglected;
                        EFTGamma2 = 0._dl
                    case (1) ! 1) Constant model;
                        if (deriv==0) then
                            EFTGamma2 = CP%EFTGamma20
                        else
                            EFTGamma2 = 0._dl
                        end if
                    case (2) ! 2) Linear model;
                        if (deriv==0) then
                            EFTGamma2= CP%EFTGamma20*a
                        else if (deriv==1) then
                            EFTGamma2= CP%EFTGamma20
                        else
                            EFTGamma2= 0._dl
                        end if
                    case (3) ! 3) Power law model;
                        if (deriv==0) then
                            EFTGamma2= CP%EFTGamma20*a**CP%EFTGamma2Exp
                        else if (deriv==1) then
                            EFTGamma2= CP%EFTGamma2Exp*CP%EFTGamma20*a**(CP%EFTGamma2Exp-1.0_dl)
                        else if (deriv==2) then
                            EFTGamma2= (CP%EFTGamma2Exp-1.0_dl)*CP%EFTGamma2Exp*CP%EFTGamma20*&
                                &a**(CP%EFTGamma2Exp-2.0_dl)
                        else if (deriv==3) then
                            EFTGamma2= (CP%EFTGamma2Exp-2.0_dl)*(CP%EFTGamma2Exp-1.0_dl)*&
                                &CP%EFTGamma2Exp*CP%EFTGamma20*a**(CP%EFTGamma2Exp-3.0_dl)
                        end if
                    case (4) ! 4) Exponential model;
                        if (deriv==0) then
                            EFTGamma2= Exp(CP%EFTGamma20*a**CP%EFTGamma2Exp)-1._dl
                        else if (deriv==1) then
                            EFTGamma2= CP%EFTGamma2Exp*CP%EFTGamma20*a**(CP%EFTGamma2Exp-1.0_dl)&
                                &*Exp(CP%EFTGamma20*a**CP%EFTGamma2Exp)
                        else if (deriv==2) then
                            EFTGamma2= a**(-2._dl+CP%EFTGamma2Exp)*Exp(a**CP%EFTGamma2Exp*CP%EFTGamma20)*&
                                &CP%EFTGamma20*CP%EFTGamma2Exp*(-1._dl+CP%EFTGamma2Exp&
                                &+a**CP%EFTGamma2Exp*CP%EFTGamma20*CP%EFTGamma2Exp)
                        else if (deriv==3) then
                            EFTGamma2= a**(-3._dl+CP%EFTGamma2Exp)*Exp(a**CP%EFTGamma2Exp*CP%EFTGamma20)*&
                                &CP%EFTGamma20*CP%EFTGamma2Exp*(2._dl-3._dl*(1._dl+a**CP%EFTGamma2Exp*&
                                &CP%EFTGamma20)*CP%EFTGamma2Exp+(1._dl+3._dl*a**CP%EFTGamma2Exp*&
                                &CP%EFTGamma20+a**(2._dl*CP%EFTGamma2Exp)*CP%EFTGamma20**2)&
                                &*CP%EFTGamma2Exp**2)
                        end if
                    case (5) ! 5) User defined;
                        stop 'You have to define your model for EFTGamma2 first! Go to the file EFT_functions.f90.'
                        if (deriv==0) then
                        else if (deriv==1) then
                        else if (deriv==2) then
                        else if (deriv==3) then
                        end if
                    case default
                        stop 'Wrong choice for EFTGamma2.'
                end select
            case (2) ! Designer mapping EFT
                select case (CP%DesignerEFTModel)
                    case (1) ! Designer f(R) models
                        EFTGamma2=0._dl
                    case (2) ! Designer mc quintessence
                        EFTGamma2=0._dl
                    case default
                        stop 'Wrong choice for EFTGamma2.'
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH Horndeski
                        EFTGamma2 = RPH_EFT_Gamma2(a, deriv)
                    case default
                        stop 'Wrong choice for EFTGamma2.'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        EFTGamma2 = 0._dl
                    case default
                        stop 'Wrong model for EFTGamma2.'
                end select
            case default
                stop 'Wrong choice for EFTGamma2.'
        end select

    end function EFTGamma2

    ! ----------------------------

    function EFTGamma3(a,deriv)
        ! EFT function Gamma_3
        implicit none
        real(dl) :: EFTGamma3
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%EFTflag)
            case (0) ! GR
                EFTGamma3 = 0._dl
            case (1) ! Pure EFT
                select case (CP%PureEFTmodelGamma3)
                    case (0) ! 0) Operator neglected;
                        EFTGamma3 = 0._dl
                    case (1) ! 1) Constant model;
                        if (deriv==0) then
                            EFTGamma3 = CP%EFTGamma30
                        else
                            EFTGamma3 = 0._dl
                        end if
                    case (2) ! 2) Linear model;
                        if (deriv==0) then
                            EFTGamma3= CP%EFTGamma30*a
                        else if (deriv==1) then
                            EFTGamma3= CP%EFTGamma30
                        else
                            EFTGamma3= 0._dl
                        end if
                    case (3) ! 3) Power law model;
                        if (deriv==0) then
                            EFTGamma3= CP%EFTGamma30*a**CP%EFTGamma3Exp
                        else if (deriv==1) then
                            EFTGamma3= CP%EFTGamma3Exp*CP%EFTGamma30*a**(CP%EFTGamma3Exp-1.0_dl)
                        else if (deriv==2) then
                            EFTGamma3= (CP%EFTGamma3Exp-1.0_dl)*CP%EFTGamma3Exp*CP%EFTGamma30*&
                                &a**(CP%EFTGamma3Exp-2.0_dl)
                        else if (deriv==3) then
                            EFTGamma3= (CP%EFTGamma3Exp-2.0_dl)*(CP%EFTGamma3Exp-1.0_dl)*&
                                &CP%EFTGamma3Exp*CP%EFTGamma30*a**(CP%EFTGamma3Exp-3.0_dl)
                        end if
                    case (4) ! 4) Exponential model;
                        if (deriv==0) then
                            EFTGamma3= Exp(CP%EFTGamma30*a**CP%EFTGamma3Exp)-1._dl
                        else if (deriv==1) then
                            EFTGamma3= CP%EFTGamma3Exp*CP%EFTGamma30*a**(CP%EFTGamma3Exp-1.0_dl)&
                                &*Exp(CP%EFTGamma30*a**CP%EFTGamma3Exp)
                        else if (deriv==2) then
                            EFTGamma3= a**(-2._dl+CP%EFTGamma3Exp)*Exp(a**CP%EFTGamma3Exp*CP%EFTGamma30)*&
                                &CP%EFTGamma30*CP%EFTGamma3Exp*(-1._dl+CP%EFTGamma3Exp&
                                &+a**CP%EFTGamma3Exp*CP%EFTGamma30*CP%EFTGamma3Exp)
                        else if (deriv==3) then
                            EFTGamma3= a**(-3._dl+CP%EFTGamma3Exp)*Exp(a**CP%EFTGamma3Exp*CP%EFTGamma30)*&
                                &CP%EFTGamma30*CP%EFTGamma3Exp*(2._dl-3._dl*(1._dl+a**CP%EFTGamma3Exp*&
                                &CP%EFTGamma30)*CP%EFTGamma3Exp+(1._dl+3._dl*a**CP%EFTGamma3Exp*&
                                &CP%EFTGamma30+a**(2._dl*CP%EFTGamma3Exp)*CP%EFTGamma30**2)&
                                &*CP%EFTGamma3Exp**2)
                        end if
                    case (5) ! 5) User defined;
                        stop 'You have to define your model for EFTGamma3 first! Go to the file EFT_functions.f90.'
                        if (deriv==0) then
                        else if (deriv==1) then
                        else if (deriv==2) then
                        else if (deriv==3) then
                        end if
                    case default
                        stop 'Wrong choice for EFTGamma3.'
                end select
            case (2) ! Designer mapping EFT
                select case (CP%DesignerEFTModel)
                    case (1) ! Designer f(R) models
                        EFTGamma3=0._dl
                    case (2) ! Designer mc quintessence
                        EFTGamma3=0._dl
                    case default
                        stop 'Wrong choice for EFTGamma3.'
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH Horndeski
                        EFTGamma3 = RPH_EFT_Gamma3(a, deriv)
                    case default
                        stop 'Wrong choice for EFTGamma3.'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        if (deriv==0) then
                            EFTGamma3 = 2._dl*(CP%Horava_lambda - CP%Horava_xi)/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                        else
                            EFTGamma3 = 0._dl
                        end if
                    case default
                        stop 'Wrong model for EFTGamma3.'
                end select
            case default
                stop 'Wrong choice for EFTGamma3.'
        end select

    end function EFTGamma3

    ! ----------------------------

    function EFTGamma4(a,deriv)
        ! EFT function Gamma_4
        implicit none
        real(dl) :: EFTGamma4
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%EFTflag)
            case (0) ! GR
                EFTGamma4 = 0._dl
            case (1) ! Pure EFT
                if (.not. CP%PureEFTHorndeski) then
                    select case (CP%PureEFTmodelGamma4)
                        case (0) ! 0) Operator neglected;
                            EFTGamma4 = 0._dl
                        case (1) ! 1) Constant model;
                            if (deriv==0) then
                                EFTGamma4 = CP%EFTGamma40
                            else
                                EFTGamma4 = 0._dl
                            end if
                        case (2) ! 2) Linear model;
                            if (deriv==0) then
                                EFTGamma4= CP%EFTGamma40*a
                            else if (deriv==1) then
                                EFTGamma4= CP%EFTGamma40
                            else
                                EFTGamma4= 0._dl
                            end if
                        case (3) ! 3) Power law model;
                            if (deriv==0) then
                                EFTGamma4= CP%EFTGamma40*a**CP%EFTGamma4Exp
                            else if (deriv==1) then
                                EFTGamma4= CP%EFTGamma4Exp*CP%EFTGamma40*a**(CP%EFTGamma4Exp-1.0_dl)
                            else if (deriv==2) then
                                EFTGamma4= (CP%EFTGamma4Exp-1.0_dl)*CP%EFTGamma4Exp*CP%EFTGamma40*&
                                    &a**(CP%EFTGamma4Exp-2.0_dl)
                            else if (deriv==3) then
                                EFTGamma4= (CP%EFTGamma4Exp-2.0_dl)*(CP%EFTGamma4Exp-1.0_dl)*&
                                    &CP%EFTGamma4Exp*CP%EFTGamma40*a**(CP%EFTGamma4Exp-3.0_dl)
                            end if
                        case (4) ! 4) Exponential model;
                            if (deriv==0) then
                                EFTGamma4= Exp(CP%EFTGamma40*a**CP%EFTGamma4Exp)-1._dl
                            else if (deriv==1) then
                                EFTGamma4= CP%EFTGamma4Exp*CP%EFTGamma40*a**(CP%EFTGamma4Exp-1.0_dl)&
                                    &*Exp(CP%EFTGamma40*a**CP%EFTGamma4Exp)
                            else if (deriv==2) then
                                EFTGamma4= a**(-2._dl+CP%EFTGamma4Exp)*Exp(a**CP%EFTGamma4Exp*CP%EFTGamma40)*&
                                    &CP%EFTGamma40*CP%EFTGamma4Exp*(-1._dl+CP%EFTGamma4Exp&
                                    &+a**CP%EFTGamma4Exp*CP%EFTGamma40*CP%EFTGamma4Exp)
                            else if (deriv==3) then
                                EFTGamma4= a**(-3._dl+CP%EFTGamma4Exp)*Exp(a**CP%EFTGamma4Exp*CP%EFTGamma40)*&
                                    &CP%EFTGamma40*CP%EFTGamma4Exp*(2._dl-3._dl*(1._dl+a**CP%EFTGamma4Exp*&
                                    &CP%EFTGamma40)*CP%EFTGamma4Exp+(1._dl+3._dl*a**CP%EFTGamma4Exp*&
                                    &CP%EFTGamma40+a**(2._dl*CP%EFTGamma4Exp)*CP%EFTGamma40**2)&
                                    &*CP%EFTGamma4Exp**2)
                            end if
                        case (5) ! 5) User defined;
                            stop 'You have to define your model for EFTGamma4 first! Go to the file EFT_functions.f90.'
                            if (deriv==0) then
                            else if (deriv==1) then
                            else if (deriv==2) then
                            else if (deriv==3) then
                            end if
                        case default
                            stop 'Wrong choice for EFTGamma4.'
                    end select
                else if (CP%PureEFTHorndeski) then
                    EFTGamma4 = -EFTGamma3(a, deriv)
                end if
            case (2) ! Designer mapping EFT
                select case (CP%DesignerEFTModel)
                    case (1) ! Designer f(R) models
                        EFTGamma4=0._dl
                    case (2) ! Designer mc quintessence
                        EFTGamma4=0._dl
                    case default
                        stop 'Wrong choice for EFTGamma4.'
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH Horndeski
                        EFTGamma4 = -RPH_EFT_Gamma3(a, deriv)
                    case default
                        stop 'Wrong choice for EFTGamma4.'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        if (deriv==0) then
                            EFTGamma4 = 2._dl*CP%Horava_xi/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                        else
                            EFTGamma4 = 0._dl
                        end if
                    case default
                        stop 'Wrong model for EFTGamma4.'
                end select
            case default
                stop 'Wrong choice for EFTGamma4.'
        end select

    end function EFTGamma4

    ! ----------------------------

    function EFTGamma5(a,deriv)
        ! EFT function Gamma_5
        implicit none
        real(dl) :: EFTGamma5
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%EFTflag)
            case (0) ! GR
                EFTGamma5 = 0._dl
            case (1) ! Pure EFT
                if (.not. CP%PureEFTHorndeski) then
                    select case (CP%PureEFTmodelGamma5)
                        case (0) ! 0) Operator neglected;
                            EFTGamma5 = 0._dl
                        case (1) ! 1) Constant model;
                            if (deriv==0) then
                                EFTGamma5 = CP%EFTGamma50
                            else
                                EFTGamma5 = 0._dl
                            end if
                        case (2) ! 2) Linear model;
                            if (deriv==0) then
                                EFTGamma5= CP%EFTGamma50*a
                            else if (deriv==1) then
                                EFTGamma5= CP%EFTGamma50
                            else
                                EFTGamma5= 0._dl
                            end if
                        case (3) ! 3) Power law model;
                            if (deriv==0) then
                                EFTGamma5= CP%EFTGamma50*a**CP%EFTGamma5Exp
                            else if (deriv==1) then
                                EFTGamma5= CP%EFTGamma5Exp*CP%EFTGamma50*a**(CP%EFTGamma5Exp-1.0_dl)
                            else if (deriv==2) then
                                EFTGamma5= (CP%EFTGamma5Exp-1.0_dl)*CP%EFTGamma5Exp*CP%EFTGamma50*&
                                    &a**(CP%EFTGamma5Exp-2.0_dl)
                            else if (deriv==3) then
                                EFTGamma5= (CP%EFTGamma5Exp-2.0_dl)*(CP%EFTGamma5Exp-1.0_dl)*&
                                    &CP%EFTGamma5Exp*CP%EFTGamma50*a**(CP%EFTGamma5Exp-3.0_dl)
                            end if
                        case (4) ! 4) Exponential model;
                            if (deriv==0) then
                                EFTGamma5= Exp(CP%EFTGamma50*a**CP%EFTGamma5Exp)-1._dl
                            else if (deriv==1) then
                                EFTGamma5= CP%EFTGamma5Exp*CP%EFTGamma50*a**(CP%EFTGamma5Exp-1.0_dl)&
                                    &*Exp(CP%EFTGamma50*a**CP%EFTGamma5Exp)
                            else if (deriv==2) then
                                EFTGamma5= a**(-2._dl+CP%EFTGamma5Exp)*Exp(a**CP%EFTGamma5Exp*CP%EFTGamma50)*&
                                    &CP%EFTGamma50*CP%EFTGamma5Exp*(-1._dl+CP%EFTGamma5Exp&
                                    &+a**CP%EFTGamma5Exp*CP%EFTGamma50*CP%EFTGamma5Exp)
                            else if (deriv==3) then
                                EFTGamma5= a**(-3._dl+CP%EFTGamma5Exp)*Exp(a**CP%EFTGamma5Exp*CP%EFTGamma50)*&
                                    &CP%EFTGamma50*CP%EFTGamma5Exp*(2._dl-3._dl*(1._dl+a**CP%EFTGamma5Exp*&
                                    &CP%EFTGamma50)*CP%EFTGamma5Exp+(1._dl+3._dl*a**CP%EFTGamma5Exp*&
                                    &CP%EFTGamma50+a**(2._dl*CP%EFTGamma5Exp)*CP%EFTGamma50**2)&
                                    &*CP%EFTGamma5Exp**2)
                            end if
                        case (5) ! 5) User defined;
                            stop 'You have to define your model for EFTGamma5 first! Go to the file EFT_functions.f90.'
                            if (deriv==0) then
                            else if (deriv==1) then
                            else if (deriv==2) then
                            else if (deriv==3) then
                            end if
                        case default
                            stop 'Wrong choice for EFTGamma5.'
                    end select
                else if (CP%PureEFTHorndeski) then
                    EFTGamma5 = +0.5_dl*EFTGamma3(a, deriv)
                end if
            case (2) ! Designer mapping EFT
                select case (CP%DesignerEFTModel)
                    case (1) ! Designer f(R) models
                        EFTGamma5=0._dl
                    case (2) ! Designer mc quintessence
                        EFTGamma5=0._dl
                    case default
                        stop 'Wrong choice for EFTGamma5.'
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH Horndeski
                        EFTGamma5 = +0.5_dl*RPH_EFT_Gamma3(a, deriv)
                    case default
                        stop 'Wrong choice for EFTGamma5.'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        EFTGamma5 = 0._dl
                    case default
                        stop 'Wrong model for EFTGamma5.'
                end select
            case default
                stop 'Wrong choice for EFTGamma5.'
        end select

    end function EFTGamma5

    ! ----------------------------

    function EFTGamma6(a,deriv)
        ! EFT function Gamma_6
        implicit none
        real(dl) :: EFTGamma6
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%EFTflag)
            case (0) ! GR
                EFTGamma6 = 0._dl
            case (1) ! Pure EFT
                if (.not. CP%PureEFTHorndeski) then
                    select case (CP%PureEFTmodelGamma6)
                        case (0) ! 0) Operator neglected;
                            EFTGamma6 = 0._dl
                        case (1) ! 1) Constant model;
                            if (deriv==0) then
                                EFTGamma6 = CP%EFTGamma60
                            else
                                EFTGamma6 = 0._dl
                            end if
                        case (2) ! 2) Linear model;
                            if (deriv==0) then
                                EFTGamma6= CP%EFTGamma60*a
                            else if (deriv==1) then
                                EFTGamma6= CP%EFTGamma60
                            else
                                EFTGamma6= 0._dl
                            end if
                        case (3) ! 3) Power law model;
                            if (deriv==0) then
                                EFTGamma6= CP%EFTGamma60*a**CP%EFTGamma6Exp
                            else if (deriv==1) then
                                EFTGamma6= CP%EFTGamma6Exp*CP%EFTGamma60*a**(CP%EFTGamma6Exp-1.0_dl)
                            else if (deriv==2) then
                                EFTGamma6= (CP%EFTGamma6Exp-1.0_dl)*CP%EFTGamma6Exp*CP%EFTGamma60*&
                                    &a**(CP%EFTGamma6Exp-2.0_dl)
                            else if (deriv==3) then
                                EFTGamma6= (CP%EFTGamma6Exp-2.0_dl)*(CP%EFTGamma6Exp-1.0_dl)*&
                                    &CP%EFTGamma6Exp*CP%EFTGamma60*a**(CP%EFTGamma6Exp-3.0_dl)
                            end if
                        case (4) ! 4) Exponential model;
                            if (deriv==0) then
                                EFTGamma6= Exp(CP%EFTGamma60*a**CP%EFTGamma6Exp)-1._dl
                            else if (deriv==1) then
                                EFTGamma6= CP%EFTGamma6Exp*CP%EFTGamma60*a**(CP%EFTGamma6Exp-1.0_dl)&
                                    &*Exp(CP%EFTGamma60*a**CP%EFTGamma6Exp)
                            else if (deriv==2) then
                                EFTGamma6= a**(-2._dl+CP%EFTGamma6Exp)*Exp(a**CP%EFTGamma6Exp*CP%EFTGamma60)*&
                                    &CP%EFTGamma60*CP%EFTGamma6Exp*(-1._dl+CP%EFTGamma6Exp&
                                    &+a**CP%EFTGamma6Exp*CP%EFTGamma60*CP%EFTGamma6Exp)
                            else if (deriv==3) then
                                EFTGamma6= a**(-3._dl+CP%EFTGamma6Exp)*Exp(a**CP%EFTGamma6Exp*CP%EFTGamma60)*&
                                    &CP%EFTGamma60*CP%EFTGamma6Exp*(2._dl-3._dl*(1._dl+a**CP%EFTGamma6Exp*&
                                    &CP%EFTGamma60)*CP%EFTGamma6Exp+(1._dl+3._dl*a**CP%EFTGamma6Exp*&
                                    &CP%EFTGamma60+a**(2._dl*CP%EFTGamma6Exp)*CP%EFTGamma60**2)&
                                    &*CP%EFTGamma6Exp**2)
                            end if
                        case (5) ! 5) User defined;
                            stop 'You have to define your model for EFTGamma6 first! Go to the file EFT_functions.f90.'
                            if (deriv==0) then
                            else if (deriv==1) then
                            else if (deriv==2) then
                            else if (deriv==3) then
                            end if
                        case default
                            stop 'Wrong choice for EFTGamma6.'
                    end select
                else if (CP%PureEFTHorndeski) then
                    EFTGamma6 = 0._dl
                end if
            case (2) ! Designer mapping EFT
                select case (CP%DesignerEFTModel)
                    case (1) ! Designer f(R) models
                        EFTGamma6=0._dl
                    case (2) ! Designer mc quintessence
                        EFTGamma6=0._dl
                    case default
                        stop 'Wrong choice for EFTGamma6.'
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH Horndeski
                        EFTGamma6 = 0._dl
                    case default
                        stop 'Wrong choice for EFTGamma6.'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        if (deriv==0) then
                            EFTGamma6 = CP%Horava_eta/4._dl/( 2._dl +2._dl*CP%Horava_xi -CP%Horava_eta )
                        else
                            EFTGamma6 = 0._dl
                        end if
                    case default
                        stop 'Wrong model for EFTGamma6.'
                end select
            case default
                stop 'Wrong choice for EFTGamma6.'
        end select

    end function EFTGamma6

    ! ----------------------------

    function EFTLambdaTemp(a,deriv)
        ! EFT function Lambda
        ! This function is defined as Lambda/Mpl*a^2
        implicit none
        real(dl):: EFTLambdaTemp
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%EFTflag)
            case (0) ! GR
                EFTLambdaTemp = 0._dl
            case (1) ! Pure EFT: Lambda is obtained through EFT designer approach (in equations.f90).
                EFTLambdaTemp = 0._dl
            case (2) ! Designer mapping EFT
                select case (CP%DesignerEFTModel)
                    case (1) ! Designer f(R) models
                        if (deriv==0) then
                            EFTLambdaTemp = Designed_EFT_Function(EFT_Lambda_Des, a, 0._dl)
                        else if (deriv==1) then
                            EFTLambdaTemp = Designed_EFT_Function(EFT_Lambdadot_Des, a, 0._dl)
                        end if
                    case (2) ! Designer mc quintessence: Lambda is obtained through EFT designer approach (in equations.f90).
                        EFTLambdaTemp = 0._dl
                    case default
                        stop 'Wrong choice for Lambda.'
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH: Lambda is obtained through EFT designer approach (in equations.f90).
                        EFTLambdaTemp = 0._dl
                    case default
                        stop 'Wrong choice for Lambda.'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        EFTLambdaTemp = 0._dl
                end select
            case default
                stop 'Wrong choice for Lambda.'
        end select

    end function EFTLambdaTemp

    ! ----------------------------

    function EFTcTemp(a,deriv)
        ! EFT function c
        ! This function is defined as c/Mpl*a^2
        implicit none
        real(dl):: EFTcTemp
        real(dl) , intent(in) :: a      ! scale factor
        integer  , intent(in) :: deriv  ! order of derivative wrt the scale factor

        select case (CP%EFTflag)
            case (0) ! GR
                EFTcTemp = 0._dl
            case (1) ! Pure EFT: c is obtained through EFT designer approach (in equations.f90).
                EFTcTemp = 0._dl
            case (2) ! Designer mapping EFT
                select case (CP%DesignerEFTModel)
                    case (1) ! Designer f(R) models
                        EFTcTemp = 0._dl
                    case (2) ! Designer mc quintessence: c is obtained through EFT designer approach (in equations.f90).
                        EFTcTemp = 0._dl
                    case default
                        stop 'Wrong choice for c.'
                end select
            case (3) ! EFT alternative parametrizations
                select case (CP%AltParEFTmodel)
                    case (1) ! RPH: c is obtained through EFT designer approach (in equations.f90).
                        EFTcTemp = 0._dl
                    case default
                        stop 'Wrong choice for c.'
                end select
            case (4) ! Full mapping EFT
                select case (CP%FullMappingEFTmodel)
                    case (1) ! Horava gravity
                        EFTcTemp = 0._dl
                end select
            case default
                stop 'Wrong choice for c.'
        end select

    end function EFTcTemp

end module EFTfunctions

! -------------------------------------------------------------------------------------------------

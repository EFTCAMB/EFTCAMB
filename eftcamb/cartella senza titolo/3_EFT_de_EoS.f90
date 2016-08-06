
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

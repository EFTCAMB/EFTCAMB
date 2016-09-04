! -------------------------------------------------------------------------------------------------
!
!   EFT_FullMapping module
!
!   Developed and implemented by:
!       Bin Hu (hu@lorentz.leidenuniv.nl)
!       Marco Raveri (mraveri@sissa.it)
!
!   EFT_FullMapping module provides EFTCAMB with consistent background (EFTw(a)) and 
!   perturbation (EFT functions, e.g. EFTOmega, EFTLambda, ...) quantities.
!
!   Input of this module are the EFTw(a) and EFT functions; 
!   The function of this module is to compute the derivatives of EFTw(a) and EFT functions;
!   Output of this module are EFTw(a) and EFT functions and their derivatives;
!
!   New developers should check the consistency of between background and perturbation quantities.
!
! -------------------------------------------------------------------------------------------------

module EFT_FullMapping
    use EFT_FullMapping_common
    use EFT_FullMapping_HS
    !use EFTDef
    implicit none

    ! 1) Definitions of the variables common to all the fullmapping code.

    real(dl), save, dimension(fm_nstep+1) :: fm_xp_pre !Sampling array before doing numerical derivatives
    real(dl), save, dimension(fm_nstep+1) :: fm_eqstate_pre,fm_wdep_pre,fm_wdepp_pre,fm_rhode_pre
    real(dl), save, dimension(fm_nstep+1) :: fm_Omega_pre,fm_Omegap_pre,fm_Omegapp_pre,fm_Omega3p_pre
    real(dl), save, dimension(fm_nstep+1) :: fm_c_pre,fm_cdot_pre
    real(dl), save, dimension(fm_nstep+1) :: fm_Lambda_pre,fm_Lambdadot_pre
    real(dl), save, dimension(fm_nstep+1) :: fm_Gamma1_pre, fm_Gamma1p_pre, fm_Gamma2_pre, fm_Gamma2p_pre
    real(dl), save, dimension(fm_nstep+1) :: fm_Gamma3_pre, fm_Gamma3p_pre, fm_Gamma4_pre, fm_Gamma4p_pre
    real(dl), save, dimension(fm_nstep+1) :: fm_Gamma5_pre, fm_Gamma5p_pre, fm_Gamma6_pre, fm_Gamma6p_pre
    real(dl), save, dimension(fm_nstep+1) :: fm_Gamma4pp_pre

    real(dl), save, dimension(fm_nstep+1-2*fm_safe) :: fm_xp   ! Array with the time sampling
    real(dl), save, dimension(fm_nstep+1-2*fm_safe) :: fm_eqstate,fm_wdep,fm_wdepp,fm_rhode
    real(dl), save, dimension(fm_nstep+1-2*fm_safe) :: fm_Omega, fm_Omegap,fm_Omegapp,fm_Omega3p
    real(dl), save, dimension(fm_nstep+1-2*fm_safe) :: fm_Lambda, fm_Lambdadot, fm_c, fm_cdot
    real(dl), save, dimension(fm_nstep+1-2*fm_safe) :: fm_Gamma1, fm_Gamma1p, fm_Gamma2, fm_Gamma2p
    real(dl), save, dimension(fm_nstep+1-2*fm_safe) :: fm_Gamma3, fm_Gamma3p, fm_Gamma4, fm_Gamma4p
    real(dl), save, dimension(fm_nstep+1-2*fm_safe) :: fm_Gamma5, fm_Gamma5p, fm_Gamma6, fm_Gamma6p
    real(dl), save, dimension(fm_nstep+1-2*fm_safe) :: fm_Gamma4pp

    real(dl), parameter :: fm_EFTturnonpiInitial = 9.d-3

contains

    ! ----------------------------------------------- !
    subroutine EFT_FM_main(success)

        implicit none

        logical, intent(inout) :: success
        integer :: fm_i
        real(dl) :: fm_err,a,h_try
        integer :: ind_deriv

        if (Feedbacklevel>1) write(*,*) 'EFTCAMB: calling full mapping code'

        if (CP%FullMappingEFTmodel .eq. 1) then
            call FM_HS_main(success) !Hu-Sawicki f(R) model
        else
            stop 'does not support for this model!'
        endif

        fm_xp_pre(:) = fm_hs_efold(:)
        fm_eqstate_pre(:) = fm_hs_eqstate(:)
        fm_rhode_pre(:) = fm_hs_rhoeff(:)
        fm_rhode_pre(:) = fm_rhode_pre(:)/3./fm_hs_m2/(1.-CP%Omegab-CP%Omegac)*(CP%Omegab+CP%Omegac)
        fm_Omega_pre(:) = fm_hs_Omega(:)
        fm_c_pre(:) = fm_hs_c(:)
        fm_Lambda_pre(:) = fm_hs_Lambda(:)

        fm_xp(:) = 0._dl
        fm_eqstate(:) = 0._dl
        fm_wdep(:) = 0._dl
        fm_wdepp(:) = 0._dl
        fm_rhode(:) = 0._dl
        fm_Omega(:) = 0._dl
        fm_Omegap(:) = 0._dl
        fm_Omegapp(:) = 0._dl
        fm_Omega3p(:) = 0._dl
        fm_Lambda(:) = 0._dl
        fm_Lambdadot(:) = 0._dl
        fm_c(:) = 0._dl
        fm_cdot(:) = 0._dl
        fm_Gamma1(:) = 0._dl
        fm_Gamma1p(:) = 0._dl
        fm_Gamma2(:) = 0._dl
        fm_Gamma2p(:) = 0._dl
        fm_Gamma3(:) = 0._dl
        fm_Gamma3p(:) = 0._dl
        fm_Gamma4(:) = 0._dl
        fm_Gamma4p(:) = 0._dl
        fm_Gamma4pp(:) = 0._dl
        fm_Gamma5(:) = 0._dl
        fm_Gamma5p(:) = 0._dl
        fm_Gamma6(:) = 0._dl
        fm_Gamma6p(:) = 0._dl

        !... prepare the wdep sampling list
        do fm_i = 1, fm_nstep+1
            a = EXP(fm_xp_pre(fm_i))
            if (a .lt. fm_EFTturnonpiInitial) then
                fm_wdep_pre(fm_i) = 0._dl
            else
                if (fm_i .le. fm_nstep) then
                    h_try = Exp(fm_xp_pre(fm_i+1)) - EXP(fm_xp_pre(fm_i))
                else
                    h_try = Exp(fm_xp_pre(fm_nstep+1))-Exp(fm_xp_pre(fm_nstep))
                endif
                h_try = h_try/2._dl
                !h_try = 0._dl
                ind_deriv = 0
                fm_wdep_pre(fm_i) = EFT_FM_dfridr(EFT_FM_WDE_Func,a,ind_deriv,h_try,fm_err)
            endif
        end do

        !... prepare the wdepp sampling list
        do fm_i = 1, fm_nstep+1
            a = EXP(fm_xp_pre(fm_i))
            if (a .lt. fm_EFTturnonpiInitial) then
                fm_wdepp_pre(fm_i) = 0._dl
            else
                if (fm_i .le. fm_nstep) then
                    h_try = Exp(fm_xp_pre(fm_i+1)) - EXP(fm_xp_pre(fm_i))
                else
                    h_try = Exp(fm_xp_pre(fm_nstep+1))-Exp(fm_xp_pre(fm_nstep))
                endif
                h_try = h_try/2._dl
                !h_try = 0._dl
                ind_deriv = 1
                fm_wdepp_pre(fm_i) = EFT_FM_dfridr(EFT_FM_WDE_Func,a,ind_deriv,h_try,fm_err)
            endif
        end do

        !... prepare the Omegap sampling list
        do fm_i = 1, fm_nstep+1
            a = EXP(fm_xp_pre(fm_i))
            if (a .lt. fm_EFTturnonpiInitial) then
                fm_Omegap_pre(fm_i) = 0._dl
            else
                if (fm_i .le. fm_nstep) then
                    h_try = Exp(fm_xp_pre(fm_i+1)) - EXP(fm_xp_pre(fm_i))
                else
                    h_try = Exp(fm_xp_pre(fm_nstep+1))-Exp(fm_xp_pre(fm_nstep))
                endif
                h_try = h_try/2._dl
                !h_try = 0._dl
                ind_deriv = 0
                fm_Omegap_pre(fm_i) = EFT_FM_dfridr(EFT_FM_Omega_Func,a,ind_deriv,h_try,fm_err)
            endif
        end do

        !... prepare the Omegapp sampling list
        do fm_i = 1, fm_nstep+1
            a = EXP(fm_xp_pre(fm_i))
            if (a .lt. fm_EFTturnonpiInitial) then
                fm_Omegapp_pre(fm_i) = 0._dl
            else
                if (fm_i .le. fm_nstep) then
                    h_try = Exp(fm_xp_pre(fm_i+1)) - EXP(fm_xp_pre(fm_i))
                else
                    h_try = Exp(fm_xp_pre(fm_nstep+1))-Exp(fm_xp_pre(fm_nstep))
                endif
                h_try= h_try/2._dl
                !h_try = 0._dl
                ind_deriv = 1
                fm_Omegapp_pre(fm_i) = EFT_FM_dfridr(EFT_FM_Omega_Func,a,ind_deriv,h_try,fm_err)
            endif
        end do

        !... prepare the Omega3p sampling list
        do fm_i = 1, fm_nstep+1
            a = EXP(fm_xp_pre(fm_i))
            if (a .lt. fm_EFTturnonpiInitial) then
                fm_Omega3p_pre(fm_i) = 0._dl
            else
                if (fm_i .le. fm_nstep) then
                    h_try = Exp(fm_xp_pre(fm_i+1)) - EXP(fm_xp_pre(fm_i))
                else
                    h_try = Exp(fm_xp_pre(fm_nstep+1))-Exp(fm_xp_pre(fm_nstep))
                endif
                h_try= h_try/2._dl
                !h_try = 0._dl
                ind_deriv = 2
                fm_Omega3p_pre(fm_i) = EFT_FM_dfridr(EFT_FM_Omega_Func,a,ind_deriv,h_try,fm_err)
            endif
        end do

        !... prepare the Lambdadot sampling list
        !... here we compute is Lambdap, not Lambdadot
        !do fm_i = 1, fm_nstep+1
        do fm_i = 1, fm_nstep+1
            a = EXP(fm_xp_pre(fm_i))
            if (a .lt. fm_EFTturnonpiInitial) then
                fm_Lambdadot_pre(fm_i) = 0._dl
            else
                if (fm_i .le. fm_nstep) then
                    h_try = Exp(fm_xp_pre(fm_i+1)) - EXP(fm_xp_pre(fm_i))
                else
                    h_try = Exp(fm_xp_pre(fm_nstep+1))-Exp(fm_xp_pre(fm_nstep))
                endif
                h_try = h_try/2._dl
                !h_try = 0._dl
                ind_deriv = 0
                fm_Lambdadot_pre(fm_i) = EFT_FM_dfridr(EFT_FM_Lambda_Func,a,ind_deriv,h_try,fm_err)
            endif
        end do

        !... copy the pre sampling list to the final one
        !... cut the the first and last fm_safe pts
        fm_eqstate(:) = fm_eqstate_pre(1+fm_safe:fm_nstep+1-fm_safe)
        fm_wdep(:) = fm_wdep_pre(1+fm_safe:fm_nstep+1-fm_safe)
        fm_wdepp(:) = fm_wdepp_pre(1+fm_safe:fm_nstep+1-fm_safe)
        fm_xp(:) = fm_xp_pre (1+fm_safe:fm_nstep+1-fm_safe)
        fm_Omega(:) = fm_Omega_pre(1+fm_safe:fm_nstep+1-fm_safe)
        fm_Omegap(:) = fm_Omegap_pre(1+fm_safe:fm_nstep+1-fm_safe)
        fm_Omegapp(:) = fm_Omegapp_pre(1+fm_safe:fm_nstep+1-fm_safe)
        fm_Omega3p(:) = fm_Omega3p_pre(1+fm_safe:fm_nstep+1-fm_safe)
        fm_c(:) = 0._dl
        fm_cdot(:) = 0._dl
        fm_Lambda(:) = fm_Lambda_pre(1+fm_safe:fm_nstep+1-fm_safe)
        fm_rhode(:) = fm_rhode_pre(1+fm_safe:fm_nstep+1-fm_safe)
        fm_Lambdadot(:) = fm_Lambdadot_pre(1+fm_safe:fm_nstep+1-fm_safe)

        !... check
        !do fm_i = 1, fm_nstep+1-2*fm_safe,1
        !  write(41,'(8E15.5)') EXP(-fm_xp(fm_i))-1.,EXP(fm_xp(fm_i)),fm_eqstate(fm_i), &
        !  & fm_c(fm_i),fm_Lambda(fm_i),fm_Omega(fm_i)
        !  write(42,'(8E15.5)') EXP(-fm_xp(fm_i))-1.,EXP(fm_xp(fm_i)),fm_Lambdadot(fm_i),&
        !  & fm_Omegap(fm_i),fm_Omegapp(fm_i),fm_Omega3p(fm_i)
        !  write(43,'(8E15.5)') EXP(-fm_xp(fm_i))-1.,EXP(fm_xp(fm_i)),fm_wdep(fm_i),fm_wdepp(fm_i),fm_rhode(fm_i)
        !end do

        !do fm_i = 1, fm_nstep+1,1
        !  write(44,'(8E15.5)') EXP(-fm_xp_pre(fm_i))-1.,EXP(fm_xp_pre(fm_i)),fm_eqstate_pre(fm_i), &
        !  & fm_c_pre(fm_i),fm_Lambda_pre(fm_i),fm_Omega_pre(fm_i)
        !  write(45,'(8E15.5)') EXP(-fm_xp_pre(fm_i))-1.,EXP(fm_xp_pre(fm_i)),fm_Lambdadot_pre(fm_i),&
        !  & fm_Omegap_pre(fm_i),fm_Omegapp_pre(fm_i),fm_Omega3p_pre(fm_i)
        !  write(46,'(8E15.5)') EXP(-fm_xp_pre(fm_i))-1.,EXP(fm_xp_pre(fm_i)),fm_wdep_pre(fm_i),fm_wdepp_pre(fm_i),fm_rhode_pre(fm_i)
        !end do

        return
    end subroutine EFT_FM_main

    ! ----------------------------------------------- !
    function EFT_FM_WDE_Func(a,ind_deriv,after)
        implicit none

        real(dl) :: a, EFT_FM_WDE_Func
        integer, intent(in) :: ind_deriv
        logical, intent(in) :: after

        if (after) then
            if (ind_deriv .eq. 0) then
                EFT_FM_WDE_Func = EFT_FM_Function_Interpolate2(fm_eqstate,a,-1._dl)
            else if (ind_deriv .eq. 1) then
                EFT_FM_WDE_Func = EFT_FM_Function_Interpolate2(fm_wdep,a,0._dl)
            else if (ind_deriv .eq. 2) then
                EFT_FM_WDE_Func = EFT_FM_Function_Interpolate2(fm_wdepp,a,0._dl)
            else if (ind_deriv .eq. 3) then
                EFT_FM_WDE_Func = a*a*EFT_FM_Function_Interpolate2(fm_rhode,a,1._dl) !bh
            else
                stop 'ind_deriv in WDEFunc is out of range in FullMapping!'
            end if
        else
            if (ind_deriv .eq. 0) then
                EFT_FM_WDE_Func = EFT_FM_Function_Interpolate(fm_eqstate_pre,a,-1._dl)
            else if (ind_deriv .eq. 1) then
                EFT_FM_WDE_Func = EFT_FM_Function_Interpolate(fm_wdep_pre,a,0._dl)
            else if (ind_deriv .eq. 2) then
                EFT_FM_WDE_Func = EFT_FM_Function_Interpolate(fm_wdepp_pre,a,0._dl)
            else if (ind_deriv .eq. 3) then
                EFT_FM_WDE_Func = a*a*EFT_FM_Function_Interpolate(fm_rhode_pre,a,1._dl) !bh
            else
                stop 'ind_deriv in WDEFunc is out of range in FullMapping!'
            end if
        endif

    end function EFT_FM_WDE_Func

    ! ----------------------------------------------- !
    function EFT_FM_Omega_Func(a,ind_deriv,after)
        implicit none

        real(dl) :: a, EFT_FM_Omega_Func
        integer, intent(in) :: ind_deriv
        logical, intent(in) :: after

        if (after) then
            if (ind_deriv .eq. 0) then
                EFT_FM_Omega_Func = EFT_FM_Function_Interpolate2(fm_Omega,a,0._dl)
            else if (ind_deriv .eq. 1) then
                EFT_FM_Omega_Func = EFT_FM_Function_Interpolate2(fm_Omegap,a,0._dl)
            else if (ind_deriv .eq. 2) then
                EFT_FM_Omega_Func = EFT_FM_Function_Interpolate2(fm_Omegapp,a,0._dl)
            else if (ind_deriv .eq. 3) then
                EFT_FM_Omega_Func = EFT_FM_Function_Interpolate2(fm_Omega3p,a,0._dl)
            else
                stop 'ind_deriv in OmegaFunc is out of range in FullMapping!'
            end if
        else
            if (ind_deriv .eq. 0) then
                EFT_FM_Omega_Func = EFT_FM_Function_Interpolate(fm_Omega_pre,a,0._dl)
            else if (ind_deriv .eq. 1) then
                EFT_FM_Omega_Func = EFT_FM_Function_Interpolate(fm_Omegap_pre,a,0._dl)
            else if (ind_deriv .eq. 2) then
                EFT_FM_Omega_Func = EFT_FM_Function_Interpolate(fm_Omegapp_pre,a,0._dl)
            else if (ind_deriv .eq. 3) then
                EFT_FM_Omega_Func = EFT_FM_Function_Interpolate(fm_Omega3p_pre,a,0._dl)
            else
                stop 'ind_deriv in OmegaFunc is out of range in FullMapping!'
            end if
        endif

    end function EFT_FM_Omega_Func

    ! ----------------------------------------------- !
    function EFT_FM_Lambda_Func(a,ind_deriv,after)
        implicit none

        real(dl) :: a, EFT_FM_Lambda_Func
        integer, intent(in) :: ind_deriv
        logical, intent(in) :: after

        if (after) then
            if (ind_deriv .eq. 0) then
                EFT_FM_Lambda_Func = EFT_FM_Function_Interpolate2(fm_Lambda,a,0._dl) !bh: need to check
            else if (ind_deriv .eq. 1) then
                EFT_FM_Lambda_Func = EFT_FM_Function_Interpolate2(fm_Lambdadot,a,0._dl)
            else
                stop 'ind_deriv in LambdaFunc is out of range in FullMapping!'
            end if
        else
            if (ind_deriv .eq. 0) then
                EFT_FM_Lambda_Func = EFT_FM_Function_Interpolate(fm_Lambda_pre,a,0._dl) !bh: need to check
            else if (ind_deriv .eq. 1) then
                EFT_FM_Lambda_Func = EFT_FM_Function_Interpolate(fm_Lambdadot_pre,a,0._dl)
            else
                stop 'ind_deriv in LambdaFunc is out of range in FullMapping!'
            end if
        endif

    end function EFT_FM_Lambda_Func

    ! ----------------------------------------------- !
    function EFT_FM_Function_Interpolate(EFTFunctionTable, a, GR_value)
        ! The fullmapping code will provide a table of sampled values for the EFT functions.
        ! This function is called in the EFT_FullMapping module to interpolate those tables.
        implicit none

        real(dl), intent(in) :: a, GR_value
        real(dl) :: EFTFunctionTable(fm_nstep+1)
        real(dl) :: EFT_FM_Function_Interpolate

        real(dl) :: x, temp, dtemp
        real(dl) :: xb(fm_ninterpol),yb(fm_ninterpol)
        integer  :: i, jlo, stint

        x = log(a)

        if (x.lt.fm_xp_pre(1)) then
            temp = GR_value
        else if (x.ge.fm_xp_pre(1) .and. x.le.fm_xp_pre(fm_nstep)) then
            ! 1) find the point corresponing to the requested input.
            call EFT_FM_hunt(fm_xp_pre, fm_nstep, x, jlo)
            ! 2) construct the table that will be interpolated
            stint = fm_ninterpol/2
            if((fm_nstep-jlo-stint).lt.0) then
                ! Last fm_ninterpol points
                do i=1, fm_ninterpol
                    xb(i)=fm_xp_pre(fm_nstep-fm_ninterpol+i)
                    yb(i)=EFTFunctionTable(fm_nstep-fm_ninterpol+i)
                end do
            else if (jlo.eq.1) then
                ! First fm_ninterpol points
                do i=1, fm_ninterpol
                    xb(i)=fm_xp_pre(i)
                    yb(i)=EFTFunctionTable(i)
                end do
            else
                ! Surrounding des_ninterpol points
                do i=1, fm_ninterpol
                    xb(i)=fm_xp_pre(jlo-stint+i)
                    yb(i)=EFTFunctionTable(jlo-stint+i)
                end do
            endif
            ! 3) Call the Neville interpolator.
            call EFT_FM_Polint(fm_ninterpol, xb, yb, x, temp, dtemp)
        else
            temp = EFTFunctionTable(fm_nstep)
        end if

        EFT_FM_Function_Interpolate = temp

    end function EFT_FM_Function_Interpolate

    ! ----------------------------------------------- !
    function EFT_FM_Function_Interpolate2(EFTFunctionTable, a, GR_value)
        ! The fullmapping code will provide a table of sampled values for the EFT functions.
        ! This function is called in the EFT_FullMapping module to interpolate those tables.
        ! This table is the final one which is directy called by eftcamb
        implicit none

        real(dl), intent(in) :: a, GR_value
        real(dl) :: EFTFunctionTable(fm_nstep+1-2*fm_safe)
        real(dl) :: EFT_FM_Function_Interpolate2

        real(dl) :: x, temp, dtemp
        real(dl) :: xb(fm_ninterpol),yb(fm_ninterpol)
        integer  :: i, jlo, stint, nsteps

        x = log(a)

        if (x.lt.fm_xp(1)) then
            temp = GR_value
        else if (x.ge.fm_xp(1) .and. x.le.fm_xp(fm_nstep-2*fm_safe)) then
            ! 1) find the point corresponing to the requested input.
            call EFT_FM_hunt(fm_xp, fm_nstep-2*fm_safe, x, jlo)

            ! 2) construct the table that will be interpolated
            stint = fm_ninterpol/2
            if((fm_nstep-2*fm_safe-jlo-stint).lt.0) then
                ! Last fm_ninterpol points
                do i=1, fm_ninterpol
                    xb(i)=fm_xp(fm_nstep-2*fm_safe-fm_ninterpol+i)
                    yb(i)=EFTFunctionTable(fm_nstep-2*fm_safe-fm_ninterpol+i)
                end do
            else if (jlo.eq.1) then
                ! First fm_ninterpol points
                do i=1, fm_ninterpol
                    xb(i)=fm_xp(i)
                    yb(i)=EFTFunctionTable(i)
                end do
            else
                ! Surrounding des_ninterpol points
                do i=1, fm_ninterpol
                    xb(i)=fm_xp(jlo-stint+i)
                    yb(i)=EFTFunctionTable(jlo-stint+i)
                end do
            endif
            ! 3) Call the Neville interpolator.
            call EFT_FM_Polint(fm_ninterpol, xb, yb, x, temp, dtemp)
        else
            temp = EFTFunctionTable(fm_nstep-2*fm_safe)
        end if

        EFT_FM_Function_Interpolate2 = temp

    end function EFT_FM_Function_Interpolate2

    ! -------------------------------------------------------------------------------------------------
    !   2) Hunting algorithm.
    !      This is used to efficiently search an ordered table by means of a hunting and
    !      a bisection algorithm.

    subroutine EFT_FM_hunt(xx,n,x,jlo)
        ! xx  = the table to be searched (in);
        ! n   = the length of the table  (in);
        ! x   = the requested value      (in);
        ! jlo = the closest (from below) entry of the table (out).
        use precision
        implicit none
        integer jlo,n
        real(dl) x,xx(n)
        integer inc,jhi,jm
        logical ascnd

        ascnd=xx(n).ge.xx(1)
        if(jlo.le.0.or.jlo.gt.n)then
            jlo=0
            jhi=n+1
            goto 3
        endif
        inc=1
        if(x.ge.xx(jlo).eqv.ascnd)then
1           jhi=jlo+inc
            if(jhi.gt.n)then
                jhi=n+1
            else if(x.ge.xx(jhi).eqv.ascnd)then
                jlo=jhi
                inc=inc+inc
                goto 1
            endif
        else
            jhi=jlo
2           jlo=jhi-inc
            if(jlo.lt.1)then
                jlo=0
            else if(x.lt.xx(jlo).eqv.ascnd)then
                jhi=jlo
                inc=inc+inc
                goto 2
            endif
        endif
3       if(jhi-jlo.eq.1)then
            if(x.eq.xx(n))jlo=n-1
            if(x.eq.xx(1))jlo=1
            return
        endif
        jm=(jhi+jlo)/2
        if(x.ge.xx(jm).eqv.ascnd)then
            jlo=jm
        else
            jhi=jm
        endif
        goto 3

    end subroutine EFT_FM_hunt

    ! -------------------------------------------------------------------------------------------------
    !   4) Neville interpolator.
    !      This is used to interpolate the EFT functions once the designer/mapping
    !      code has sampled them

    subroutine EFT_FM_Polint(n,xa,ya,xpl,ypl,dypl)
        ! n    = number of points in the table                          (in)
        ! xa   = first coordinate of the points to be interpolated      (in)
        ! ya   = second coordinate of the points to be interpolated     (in)
        ! xpl  = requested value of x                                   (in)
        ! ypl  = value of the interpolated function at xpl              (out)
        ! dypl =
        use precision
        implicit none

        integer , intent(in) :: n               ! Length of the table of points
        real(dl), intent(in) :: xa(n),ya(n)     ! The tables of points to be interpolated as y=f(x)
        real(dl) :: dypl,xpl,ypl

        integer  :: i,m,ns
        real(dl) :: den,dif,dift,ho,hp,wpl,cc(n),d(n)

        ns=1
        dif=abs(xpl-xa(1))
        do i=1,n
            dift=abs(xpl-xa(i))
            if (dift<dif) then
                ns=i
                dif=dift
            endif
            cc(i)=ya(i)
            d(i)=ya(i)
        end do
        ypl=ya(ns)
        ns=ns-1
        do m=1,n-1
            do i=1,n-m
                ho=xa(i)-xpl
                hp=xa(i+m)-xpl
                wpl=cc(i+1)-d(i)
                den=ho-hp
                if (den==0.) then
                    write(*,*) 'failure in polint in FullMapping'
                    stop
                end if
                den=wpl/den
                d(i)=hp*den
                cc(i)=ho*den
            end do
            if (2*ns<n-m) then
                dypl=cc(ns+1)
            else
                dypl=d(ns)
                ns=ns-1
            endif
            ypl=ypl+dypl
        end do

        return
    end subroutine EFT_FM_Polint

    ! ----------------------------------------------- !
    ! algorithm to compute the numerical derivative.
    function EFT_FM_dfridr(func,x,ind,h,err)

        implicit none

        integer,parameter :: ntab = 100
        integer, intent(in) :: ind
        real(dl) :: EFT_FM_dfridr,err,h,x,func
        external func

        real(dl), parameter :: CON=1.4_dl       ! decrease of the stepsize.
        real(dl), parameter :: CON2=CON*CON
        real(dl), parameter :: BIG=1.d+30
        real(dl), parameter :: SAFE=2._dl

        integer i,j
        real(dl) errt, fac, hh
        real(dl), dimension(ntab,ntab) :: a

        if (h.eq.0._dl) h = 1.d-8

        hh=h
        a(1,1)=(func(x+hh,ind,.false.)-func(x-hh,ind,.false.))/(2.0_dl*hh)
        err=BIG

        EFT_FM_dfridr = 0._dl

        do 12 i=2,NTAB
            hh=hh/CON
            a(1,i)=(func(x+hh,ind,.false.)-func(x-hh,ind,.false.))/(2.0_dl*hh)
            fac=CON2
            do 11 j=2,i
                a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1._dl)
                fac=CON2*fac
                errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                if (errt.le.err) then
                    err=errt
                    EFT_FM_dfridr=a(j,i)
                endif
11          continue

            if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err) return

12      continue

        return
    end function EFT_FM_dfridr

end module EFT_FullMapping
!--------------------------------------------------------------------------------------------!



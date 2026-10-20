! Verifier driver (t10a-fullmap-hdotdotdot-formula): call the COMPILED full-mapping compute_H_derivs
! (EFTCAMBFullModelComputeHubbleDer, reached through a concrete EFTCAMB_Kmouflage instance) on input states
! read from stdin, and print Hdot, Hdotdot, Hdotdotdot, plus the finding's proposed Hdotdotdot expression
! transcribed verbatim from the finding text (so the fix is tested on the same inputs).
! Input lines: a H O0 O1 O2 O3 O4 P Pd Pdd L Ld Ldd   (P = gpresm_t, L = EFTLambda, O_k = d^k Omega/da^k)
! Output lines: Hdot Hdotdot Hdotdotdot Hdotdotdot_proposed
program t10a_hddd_toy_driver

    use precision
    use EFTCAMB_cache
    use EFTCAMB_Kmouflage_Mod

    implicit none

    type(EFTCAMB_Kmouflage) :: model
    type(TEFTCAMB_parameter_cache) :: par_cache
    type(TEFTCAMB_timestep_cache) :: eft_cache
    real(dl) :: a, H, O0, O1, O2, O3, O4, P, Pd, Pdd, L, Ld, Ldd
    real(dl) :: a2, adotoa2, one_plus_omega, denom, N1, prop
    integer :: ios

    do
        read(*, *, iostat=ios) a, H, O0, O1, O2, O3, O4, P, Pd, Pdd, L, Ld, Ldd
        if (ios /= 0) exit
        call eft_cache%initialize()
        eft_cache%a = a
        eft_cache%adotoa = H
        eft_cache%EFTOmegaV = O0
        eft_cache%EFTOmegaP = O1
        eft_cache%EFTOmegaPP = O2
        eft_cache%EFTOmegaPPP = O3
        eft_cache%EFTOmegaPPPP = O4
        eft_cache%gpresm_t = P
        eft_cache%gpresdotm_t = Pd
        eft_cache%gpresdotdotm_t = Pdd
        eft_cache%EFTLambda = L
        eft_cache%EFTLambdadot = Ld
        eft_cache%EFTLambdadotdot = Ldd

        call model%compute_H_derivs(a, par_cache, eft_cache)

        ! finding's proposed replacement, verbatim (uses the routine's own Hdot, Hdotdot)
        a2 = a*a
        adotoa2 = eft_cache%adotoa**2
        one_plus_omega = 1.0_dl + eft_cache%EFTOmegaV
        denom = one_plus_omega + 0.5_dl*a*eft_cache%EFTOmegaP
        N1 = a*eft_cache%adotoa*(1.5_dl*eft_cache%EFTOmegaP + 0.5_dl*a*eft_cache%EFTOmegaPP)
        prop = 1._dl/denom*( -N1*eft_cache%Hdotdot &
            - a*adotoa2*(0.5_dl*adotoa2 + 1.5_dl*eft_cache%Hdot)*(3._dl*eft_cache%EFTOmegaP + &
            4._dl*a*eft_cache%EFTOmegaPP + a2*eft_cache%EFTOmegaPPP) &
            - 0.5_dl*a2*eft_cache%adotoa**4*(7._dl*eft_cache%EFTOmegaPP + 6._dl*a*eft_cache%EFTOmegaPPP + &
            a2*eft_cache%EFTOmegaPPPP) &
            - (eft_cache%Hdot**2 + eft_cache%adotoa*eft_cache%Hdotdot)*(one_plus_omega + &
            3.5_dl*a*eft_cache%EFTOmegaP + 1.5_dl*a2*eft_cache%EFTOmegaPP) &
            - 0.5_dl*a*adotoa2*eft_cache%Hdot*(9._dl*eft_cache%EFTOmegaP + 13._dl*a*eft_cache%EFTOmegaPP + &
            3._dl*a2*eft_cache%EFTOmegaPPP) &
            - 0.5_dl*(eft_cache%gpresdotdotm_t + 4._dl*eft_cache%adotoa*eft_cache%gpresdotm_t + &
            2._dl*eft_cache%Hdot*eft_cache%gpresm_t + 4._dl*adotoa2*eft_cache%gpresm_t) &
            - 0.5_dl*(eft_cache%EFTLambdadotdot + 4._dl*eft_cache%adotoa*eft_cache%EFTLambdadot + &
            2._dl*eft_cache%Hdot*eft_cache%EFTLambda + 4._dl*adotoa2*eft_cache%EFTLambda))

        write(*, '(4ES26.17)') eft_cache%Hdot, eft_cache%Hdotdot, eft_cache%Hdotdotdot, prop
    end do

end program t10a_hddd_toy_driver

"""Apply the finding's proposed Hdotdotdot formula to a tree copy (06p1) and, on the python get_eft_functions path
(11), supply EFTLambdadotdot / EFTOmegaPPPP by central differences of the model's background EFT functions before the
third-derivative compute_H_derivs call (full models only). Usage: python patch_hdddot_fix.py <tree copy>"""

import sys
from pathlib import Path

tree = Path(sys.argv[1])
assert "work/verify" in str(tree), "only patch verifier copies"
f6 = tree / "fortran/eftcamb/06p1_abstract_EFTCAMB_full_map.f90"
s = f6.read_text()
i0 = s.index("    denom = one_plus_omega + 0.5*a*eft_cache%EFTOmegaP")
i1 = s.index("    end subroutine EFTCAMBFullModelComputeHubbleDer")
new = """    denom = one_plus_omega + 0.5_dl*a*eft_cache%EFTOmegaP
    denomdot = a*eft_cache%adotoa*(1.5_dl*eft_cache%EFTOmegaP + 0.5_dl*a*eft_cache%EFTOmegaPP)
    S3 = 3._dl*eft_cache%EFTOmegaP + 4._dl*a*eft_cache%EFTOmegaPP + a2*eft_cache%EFTOmegaPPP
    T3 = one_plus_omega + 3.5_dl*a*eft_cache%EFTOmegaP + 1.5_dl*a2*eft_cache%EFTOmegaPP
    N1 = -0.5_dl*a*adotoa2*(adotoa2 + 3._dl*eft_cache%Hdot)*S3 &
        - 0.5_dl*a2*adotoa2**2*(7._dl*eft_cache%EFTOmegaPP + 6._dl*a*eft_cache%EFTOmegaPPP + a2*eft_cache%EFTOmegaPPPP) &
        - (eft_cache%Hdot**2 + eft_cache%adotoa*eft_cache%Hdotdot)*T3 &
        - 0.5_dl*a*adotoa2*eft_cache%Hdot*(9._dl*eft_cache%EFTOmegaP + 13._dl*a*eft_cache%EFTOmegaPP + &
        3._dl*a2*eft_cache%EFTOmegaPPP) &
        - 0.5_dl*(eft_cache%gpresdotdotm_t + 4._dl*eft_cache%adotoa*eft_cache%gpresdotm_t + &
        2._dl*eft_cache%Hdot*eft_cache%gpresm_t + 4._dl*adotoa2*eft_cache%gpresm_t) &
        - 0.5_dl*(eft_cache%EFTLambdadotdot + 4._dl*eft_cache%adotoa*eft_cache%EFTLambdadot + &
        2._dl*eft_cache%Hdot*eft_cache%EFTLambda + 4._dl*adotoa2*eft_cache%EFTLambda)
    eft_cache%Hdotdotdot = (N1 - denomdot*eft_cache%Hdotdot)/denom

"""
s = s[:i0] + new + s[i1:]
s = s.replace("    real(dl) :: a2, denom, one_plus_omega, adotoa2, N1\n",
              "    real(dl) :: a2, denom, one_plus_omega, adotoa2, N1, denomdot, S3, T3\n", 1)
f6.write_text(s)

f11 = tree / "fortran/eftcamb/11_EFTCAMB_python.f90"
s = f11.read_text()
s = s.replace("    real(dl) :: a, adotoa, grhorad\n",
              "    real(dl) :: a, adotoa, grhorad\n    type(TEFTCAMB_timestep_cache) :: tcp, tcm\n"
              "    real(dl), parameter :: hfd = 1.e-4_dl\n", 1)
anchor = "                ! third derivative of H:\n                call EFTCAMB%model%compute_H_derivs(a, par_cache, cache)\n"
assert s.count(anchor) == 1
s = s.replace(anchor, """                ! VERIFIER PATCH: inputs of Hdotdotdot not computed on this path for full models
                if (.not. EFTCAMB%EFTCAMB_model_is_designer) then
                    tcp = cache
                    tcm = cache
                    call EFTCAMB%model%compute_background_EFT_functions(a*(1._dl + hfd), par_cache, tcp)
                    call EFTCAMB%model%compute_background_EFT_functions(a*(1._dl - hfd), par_cache, tcm)
                    cache%EFTOmegaPPPP = (tcp%EFTOmegaPPP - tcm%EFTOmegaPPP)/(2._dl*hfd*a)
                    cache%EFTLambdadotdot = a*adotoa*(tcp%EFTLambdadot - tcm%EFTLambdadot)/(2._dl*hfd*a) - &
                        2._dl*adotoa*cache%EFTLambdadot
                end if
""" + anchor, 1)
f11.write_text(s)
print("patched", f6, f11)

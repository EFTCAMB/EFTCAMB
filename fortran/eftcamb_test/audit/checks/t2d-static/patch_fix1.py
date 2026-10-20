"""Apply the candidate fix for compute_additional_derivs (designer + full mapping) to a TREE COPY.

Usage: python patch_fix1.py TREE_COPY
Inside dummy_helper of both files: evaluate everything at temp_a into temp_cache, return EFTcdot / EFTLambdadot
(so that dfridr gives d(Xdot)/da and Xdotdot = a*H*dXdot/da - 2*H*Xdot), set gamma2PP/gamma3PP/gamma3PPP from the
temporary cache, and add EFTLambdadotdot to the designer version.
"""

import re
import sys
from pathlib import Path

tree = Path(sys.argv[1])


def fix_helper(src: str, start_pat: str, end_pat: str, repl) -> str:
    i = src.index(start_pat)
    j = src.index(end_pat, i)
    return src[:i] + repl(src[i:j]) + src[j:]


def common(block: str) -> str:
    block = block.replace("call self%compute_adotoa(a, temp_par_cache, eft_cache)",
                          "call self%compute_adotoa(temp_a, temp_par_cache, temp_cache)")
    block = block.replace("call self%compute_background_EFT_functions(a, temp_par_cache, eft_cache)",
                          "call self%compute_background_EFT_functions(temp_a, temp_par_cache, temp_cache)")
    block = block.replace("temp_par_cache%grhormass(nu_i)/a**2", "temp_par_cache%grhormass(nu_i)/temp_a**2")
    block = block.replace("(a*temp_par_cache%nu_masses(nu_i)", "(temp_a*temp_par_cache%nu_masses(nu_i)")
    block = block.replace("eft_cache%Hdot", "temp_cache%Hdot").replace("temp_temp_cache", "temp_cache")
    block = block.replace("cP = eft_cache%EFTcdot/(a*eft_cache%adotoa) + 2._dl*eft_cache%EFTc/a",
                          "cP = temp_cache%EFTcdot")
    block = block.replace("LambdaP = eft_cache%EFTLambdadot/(a*eft_cache%adotoa) + 2._dl*eft_cache%EFTLambda/a",
                          "LambdaP = temp_cache%EFTLambdadot")
    return block


# designer
p = tree / "fortran/eftcamb/06p2_abstract_EFTCAMB_designer.f90"
s = p.read_text()
s = fix_helper(s, "subroutine dummy_helper(temp_a", "end subroutine dummy_helper", common)
s = s.replace(
    "    eft_cache%EFTcdotdot = a*eft_cache%adotoa*dEFTcdot_da - 2._dl*eft_cache%adotoa*eft_cache%EFTcdot\n\n    contains",
    "    eft_cache%EFTcdotdot = a*eft_cache%adotoa*dEFTcdot_da - 2._dl*eft_cache%adotoa*eft_cache%EFTcdot\n"
    "    eft_cache%EFTLambdadotdot = a*eft_cache%adotoa*dfridr(dummy_helper_LambdaP, a, 0.03_dl*a, err) - &\n"
    "        2._dl*eft_cache%adotoa*eft_cache%EFTLambdadot\n\n    contains")
s = s.replace(
    "    end subroutine EFTCAMBDesignerModelComputeAdditionalDer",
    "    function dummy_helper_LambdaP(a_temp)\n    real(dl) :: a_temp\n    real(dl) :: dummy_helper_LambdaP\n"
    "    real(dl) :: gamma1P, gamma2P, gamma3P, cP\n    call dummy_helper(a_temp, gamma1P, gamma2P, gamma3P, cP)\n"
    "    dummy_helper_LambdaP = cP\n    end function dummy_helper_LambdaP\n\n"
    "    end subroutine EFTCAMBDesignerModelComputeAdditionalDer")
# LambdaP for designer: reuse helper output slot through a module-level switch is overkill; return Lambdadot via cP
# when called from dummy_helper_LambdaP. Implement by a host logical.
s = s.replace("    real(dl) :: dEFTcdot_da, err\n", "    real(dl) :: dEFTcdot_da, err\n    logical :: want_lambda\n\n    want_lambda = .false.\n", 1)
s = s.replace("    cP = temp_cache%EFTcdot\n", "    cP = temp_cache%EFTcdot\n    if (want_lambda) cP = temp_cache%EFTLambdadot\n", 1)
s = s.replace("    call dummy_helper(a_temp, gamma1P, gamma2P, gamma3P, cP)\n    dummy_helper_LambdaP = cP",
              "    want_lambda = .true.\n    call dummy_helper(a_temp, gamma1P, gamma2P, gamma3P, cP)\n    want_lambda = .false.\n"
              "    dummy_helper_LambdaP = cP")
p.write_text(s)

# full mapping
p = tree / "fortran/eftcamb/06p1_abstract_EFTCAMB_full_map.f90"
s = p.read_text()


def full(block: str) -> str:
    block = common(block)
    block = block.replace("    gamma3P = temp_cache%EFTGamma3P\n",
                          "    gamma3P = temp_cache%EFTGamma3P\n    gamma2PP = temp_cache%EFTGamma2PP\n"
                          "    gamma3PP = temp_cache%EFTGamma3PP\n    gamma3PPP = temp_cache%EFTGamma3PPP\n")
    return block


s = fix_helper(s, "subroutine dummy_helper(temp_a", "end subroutine dummy_helper", full)
p.write_text(s)
print("patched", tree)

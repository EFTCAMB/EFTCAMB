"""Verifier build for finding minkowski-tachpos-lambdadotdot-zero: copy of A0 in which EFTStabilityComputation supplies
EFTLambdadotdot (when the model left it at 0) BEFORE compute_positivity_bounds, from a central second difference of the
model's own EFTLambda VALUES in ln a (Lambda = EFTLambda/a^2, EFTLambdadotdot = a^2 (Hdot Lambda_N + H^2 Lambda_NN)),
step h = 0.08.  Diagnostic only (not a proposed fix): shows the end-to-end EFT_minkowski_limit verdict once the
coefficient is no longer zero.

Usage: python minktach_build_patch.py
"""

import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ROOT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
            "scratchpad/audit")
DEST = ROOT / "work/verify/minktach/P_val"

if not DEST.exists():
    shutil.copytree(ROOT / "trees/A0", DEST, symlinks=True)
src = DEST / "fortran/eftcamb/09_EFTCAMB_stability.f90"
txt = src.read_text()
if "VERIFY PATCH" not in txt:
    head, sep, tail = txt.partition("    subroutine EFTStabilityComputation(")
    assert sep
    decl_old = "    integer  :: nu_i, ind, ind_max\n"
    assert decl_old in tail
    tail = tail.replace(decl_old, decl_old + "    type(TEFTCAMB_timestep_cache) :: vfy_p, vfy_m  ! VERIFY PATCH\n"
                        "    real(dl) :: vfy_h, vfy_lp, vfy_lm, vfy_l0, vfy_ln, vfy_lnn  ! VERIFY PATCH\n", 1)
    call_old = "    ! Compute positivity bounds:\n    call input_model%compute_positivity_bounds("
    assert call_old in tail
    patch = (
        "    ! VERIFY PATCH: EFTLambdadotdot from EFTLambda values when the model leaves it at zero\n"
        "    if ((.not. input_EFTCAMB%EFTCAMB_model_is_designer) .and. eft_cache%EFTLambdadotdot == 0._dl) then\n"
        "        vfy_h = 0.08_dl\n"
        "        vfy_p = eft_cache\n"
        "        vfy_m = eft_cache\n"
        "        call input_model%compute_background_EFT_functions(a*exp(vfy_h), params_cache, vfy_p)\n"
        "        call input_model%compute_background_EFT_functions(a*exp(-vfy_h), params_cache, vfy_m)\n"
        "        vfy_lp = vfy_p%EFTLambda/(a*exp(vfy_h))**2\n"
        "        vfy_lm = vfy_m%EFTLambda/(a*exp(-vfy_h))**2\n"
        "        vfy_l0 = eft_cache%EFTLambda/a**2\n"
        "        vfy_ln = (vfy_lp - vfy_lm)/(2._dl*vfy_h)\n"
        "        vfy_lnn = (vfy_lp - 2._dl*vfy_l0 + vfy_lm)/vfy_h**2\n"
        "        eft_cache%EFTLambdadotdot = a**2*(eft_cache%Hdot*vfy_ln + eft_cache%adotoa**2*vfy_lnn)\n"
        "    end if\n"
    )
    tail = tail.replace(call_old, patch + call_old, 1)
    src.write_text(head + sep + tail)
res = harness.build_tree(ROOT / "trees/A0", DEST)
print(json.dumps(res, indent=1))

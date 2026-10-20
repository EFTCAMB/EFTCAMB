"""Build a copy of A0 in which EFTStabilityComputation supplies EFTLambdadotdot for non-designer (full-mapping) models
before compute_positivity_bounds, by a central difference (in ln a) of the model's own EFTLambdadot from
compute_background_EFT_functions at a(1 +- h).  Nothing else changes.  Used only to show the end-to-end effect of the
unset EFTLambdadotdot on the EFT_minkowski_limit verdict (not a proposed fix).

Usage: python hderivs_build_ldd_patch.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ROOT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
            "scratchpad/audit")
DEST = ROOT / "work/verify/hderivs/P_ldd"

if not DEST.exists():
    import shutil

    shutil.copytree(ROOT / "trees/A0", DEST, symlinks=True)
src = DEST / "fortran/eftcamb/09_EFTCAMB_stability.f90"
txt = src.read_text()
if "AUDIT PATCH" not in txt:
    head, sep, tail = txt.partition("    subroutine EFTStabilityComputation(")
    decl_old = "    integer  :: nu_i, ind, ind_max\n"
    assert decl_old in tail
    tail = tail.replace(decl_old, decl_old + "    type(TEFTCAMB_timestep_cache) :: audit_tp, audit_tm  ! AUDIT PATCH\n"
                        "    real(dl) :: audit_h  ! AUDIT PATCH\n", 1)
    call_old = "    ! Compute positivity bounds:\n    call input_model%compute_positivity_bounds("
    assert call_old in tail
    patch = (
        "    ! AUDIT PATCH: EFTLambdadotdot from a central difference of the model EFTLambdadot (full-mapping models)\n"
        "    if (.not. input_EFTCAMB%EFTCAMB_model_is_designer) then\n"
        "        audit_h = 1.e-3_dl\n"
        "        audit_tp = eft_cache\n"
        "        audit_tm = eft_cache\n"
        "        call input_model%compute_background_EFT_functions(a*(1._dl + audit_h), params_cache, audit_tp)\n"
        "        call input_model%compute_background_EFT_functions(a*(1._dl - audit_h), params_cache, audit_tm)\n"
        "        eft_cache%EFTLambdadotdot = eft_cache%adotoa*(audit_tp%EFTLambdadot - audit_tm%EFTLambdadot)/ &\n"
        "            (log(1._dl + audit_h) - log(1._dl - audit_h)) - 2._dl*eft_cache%adotoa*eft_cache%EFTLambdadot\n"
        "    end if\n"
    )
    tail = tail.replace(call_old, patch + call_old, 1)
    src.write_text(head + sep + tail)
    print("patched", src)
r = harness.build_tree(ROOT / "trees/A0", DEST)
print(json.dumps(r, indent=1))

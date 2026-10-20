"""Verifier (pi-rk45-step-instability): build an A0 copy whose scalar RK45 max step is capped to a constant
AUDIT_VFY_MAXSTEP Mpc (env var, read once) while EFT is active. Independent of the omega_pi-based knob."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

DEST = harness.AUDIT_ROOT / "work" / "verify" / "pistep_cap"
OLD = "        rk_settings%max_step_size = 20._dl/CP%Accuracy%IntTolBoost\n    end if\n"
NEW = OLD + """    block
        character(len=32) :: vfy_buf
        real(dl), save :: vfy_cap = -1._dl
        integer :: vfy_len, vfy_st
        if (vfy_cap < 0._dl) then
            call get_environment_variable("AUDIT_VFY_MAXSTEP", vfy_buf, vfy_len, vfy_st)
            vfy_cap = 0._dl
            if (vfy_st == 0 .and. vfy_len > 0) read(vfy_buf, *) vfy_cap
        end if
        if (vfy_cap > 0._dl .and. CP%EFTCAMB%EFTflag /= 0 .and. EV%EFTCAMBactive) &
            rk_settings%max_step_size = min(rk_settings%max_step_size, vfy_cap)
    end block
"""
if not DEST.exists():
    import shutil
    shutil.copytree(harness.tree_path("A0"), DEST, symlinks=True)
f = DEST / "fortran" / "equations.f90"
s = f.read_text()
if "AUDIT_VFY_MAXSTEP" not in s:
    assert s.count(OLD) == 1
    f.write_text(s.replace(OLD, NEW))
print(harness.build_tree("A0", DEST))

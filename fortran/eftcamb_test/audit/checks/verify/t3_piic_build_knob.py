"""Verifier (t3-pi-ic-singular-denominator): build an A0 copy whose EFT_IC_type=1 pi IC is regularised by env knob
AUDIT_VPI_REG=f: pi = pi' = 0 when |C + k^2 D1 + k^4 D2| < f (|C| + k^2 |D1| + k^4 |D2|). f unset/0 -> unchanged."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

DEST = harness.AUDIT_ROOT / "work" / "verify" / "t3piic_knob"
src = harness.tree_path("A0")
if not DEST.exists():
    import shutil
    shutil.copytree(src, DEST, symlinks=True)
f = DEST / "fortran" / "eftcamb" / "09_EFTCAMB_IC.f90"
s = f.read_text()
old = "            if (EV%eft_cache%EFTpiC + k2*EV%eft_cache%EFTpiD1 + k2**2*EV%eft_cache%EFTpiD2 /= 0._dl) then\n"
new = ("            call get_environment_variable('AUDIT_VPI_REG', vpi_env, status=vpi_stat)\n"
       "            vpi_f = 0._dl\n"
       "            if (vpi_stat == 0) read(vpi_env, *) vpi_f\n"
       "            if (abs(EV%eft_cache%EFTpiC + k2*EV%eft_cache%EFTpiD1 + k2**2*EV%eft_cache%EFTpiD2) > &\n"
       "                vpi_f*(abs(EV%eft_cache%EFTpiC) + k2*abs(EV%eft_cache%EFTpiD1) + k2**2*abs(EV%eft_cache%EFTpiD2)) &\n"
       "                .and. EV%eft_cache%EFTpiC + k2*EV%eft_cache%EFTpiD1 + k2**2*EV%eft_cache%EFTpiD2 /= 0._dl) then\n")
if "AUDIT_VPI_REG" not in s:
    assert s.count(old) == 1
    s = s.replace(old, new)
    # declarations: insert after the subroutine's first 'implicit none'? find the routine containing the case
    i = s.index("AUDIT_VPI_REG")
    j = s.rfind("implicit none", 0, i)
    k = s.index("\n", j) + 1
    s = s[:k] + "    character(len=64) :: vpi_env\n    integer :: vpi_stat\n    real(dl) :: vpi_f\n" + s[k:]
    f.write_text(s)
print(harness.build_tree(DEST, DEST))

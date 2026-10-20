"""Build an instrumented release copy of A0: every equispaced-interpolation index computation reports when ind >= N
(printing x, ind, x(N), and the out-of-bounds x(N+1), y(N+1) actually read), and clamps ind to N-1 only when the
environment variable EFT_AUDIT_CLAMP=1 is set, so one binary gives both the as-shipped and the clamped behaviour."""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness  # noqa: E402

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
DEST = AUDIT / "work/verify/interp_xfinal_instr"
if not DEST.exists():
    import shutil

    shutil.copytree(AUDIT / "trees/A0", DEST, symlinks=True)
src = DEST / "fortran/eftcamb/02_equispaced_interpolation_linear_1D.f90"
txt = src.read_text()
if "audit_oob_hit" not in txt:
    pat = re.compile(r"^(\s*)ind = int\(\(x - self%x_initial\)/self%grid_width\) \+ 1\n", re.M)
    txt, n = pat.subn(lambda m: m.group(0) + f"{m.group(1)}if (ind >= self%num_points) call audit_oob_hit(self, x, ind)\n", txt)
    assert n == 7, n
    helper = """
    subroutine audit_oob_hit(self, x, ind)
    class(equispaced_linear_interpolate_function_1D) :: self
    real(dl), intent(in) :: x
    integer, intent(inout) :: ind
    integer, save :: clamp = -1, nhit = 0
    character(len=8) :: val
    integer :: n
    if (clamp < 0) then
        call get_environment_variable('EFT_AUDIT_CLAMP', val)
        clamp = 0
        if (trim(val) == '1') clamp = 1
    end if
    n = self%num_points
    nhit = nhit + 1
    if (nhit <= 400) write(0, '(a,i5,a,es24.16,a,i6,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,i1)') &
        'AUDIT_OOB hit=', nhit, ' x=', x, ' ind=', ind, ' x_final=', self%x_final, ' x(N)=', self%x(n), &
        ' x(N+1)=', self%x(n + 1), ' y(N+1)=', self%y(n + 1), ' clamp=', clamp
    if (clamp == 1) ind = n - 1
    end subroutine audit_oob_hit

    end module equispaced_linear_interpolation_1D"""
    txt, n = re.subn(r"\n\s*end module equispaced_linear_interpolation_1D", helper, txt)
    assert n == 1
    src.write_text(txt)
print(harness.build_tree(DEST, DEST))

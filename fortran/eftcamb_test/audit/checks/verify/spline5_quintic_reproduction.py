"""Independent check of finding spline5-poly5-dx6.

Compiles the repo's 02p12_splines5.f90 (copied verbatim to a scratch dir, plus a patched copy with
c(6)*dx**5) into a tiny Fortran test program. A quintic spline with bctype=[2,2] ("derivatives consistent
with interpolating quintic") must reproduce a global quintic polynomial exactly. Prints the max error of
poly5 vs the exact polynomial and the consistency of dpoly5 with a finite difference of poly5.
"""

import os
import shutil
import subprocess

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))
SRC = os.path.join(REPO, "fortran/eftcamb/02_bundled_algorithms/02p12_splines5.f90")
WORK = (
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
    "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/spline5"
)

PREC = """module precision
implicit none
integer, parameter :: dl = kind(1.d0)
end module precision
"""

PROG = """program t
use precision
use splines5
implicit none
integer, parameter :: n = 9
real(dl) :: x(n), y(n), c(0:6, n-1), xx, e, emax, emaxnode, fd, d, edmax, h
integer :: i, j
x = [0.01_dl, 0.05_dl, 0.1_dl, 0.25_dl, 0.4_dl, 0.6_dl, 0.8_dl, 0.9_dl, 1.0_dl]
do i = 1, n
  y(i) = f(x(i))
end do
call spline5pars(x, y, [2, 2], [0._dl, 0._dl, 0._dl, 0._dl], c)
emax = 0; emaxnode = 0; edmax = 0
do i = 1, n-1
  ! right end of segment i evaluated with segment i's polynomial
  e = abs(poly5(x(i+1), c(:, i)) - y(i+1))
  emaxnode = max(emaxnode, e)
  do j = 0, 20
    xx = x(i) + (x(i+1)-x(i))*j/20._dl
    emax = max(emax, abs(poly5(xx, c(:, i)) - f(xx)))
    h = 1.e-6_dl*(x(i+1)-x(i))
    fd = (poly5(xx+h, c(:, i)) - poly5(xx-h, c(:, i)))/(2*h)
    d = dpoly5(xx, c(:, i))
    edmax = max(edmax, abs(d - fd)/max(abs(fd), 1.e-30_dl))
  end do
end do
write(*, '(a,es12.4)') 'max |poly5 - exact quintic| on segments  = ', emax
write(*, '(a,es12.4)') 'max |poly5(x_{i+1}, seg i) - y_{i+1}|     = ', emaxnode
write(*, '(a,es12.4)') 'max |dpoly5 - FD(poly5)|/|FD|            = ', edmax
contains
real(dl) function f(z)
real(dl), intent(in) :: z
f = 0.3_dl - 1.1_dl*z + 2._dl*z**2 - 0.7_dl*z**3 + 1.5_dl*z**4 - 2.2_dl*z**5
end function
end program
"""


def build_and_run(tag, patch):
    d = os.path.join(WORK, tag)
    os.makedirs(d, exist_ok=True)
    src = open(SRC).read()
    if patch:
        assert src.count("c(6)*dx**6") == 1
        src = src.replace("c(6)*dx**6", "c(6)*dx**5")
    open(os.path.join(d, "splines5.f90"), "w").write(src)
    open(os.path.join(d, "precision.f90"), "w").write(PREC)
    open(os.path.join(d, "t.f90"), "w").write(PROG)
    gf = shutil.which("gfortran") or "/Users/marco/Software/miniforge3/bin/gfortran"
    subprocess.run(
        [gf, "-O0", "-g", "-fcheck=all", "-ffree-line-length-none", "precision.f90", "splines5.f90", "t.f90",
         "-llapack", "-lblas", "-o", "t"],
        cwd=d, check=True, env={**os.environ, "LIBRARY_PATH": "/Users/marco/Software/miniforge3/lib"},
    )
    out = subprocess.run(["./t"], cwd=d, capture_output=True, text=True,
                         env={**os.environ, "DYLD_LIBRARY_PATH": "/Users/marco/Software/miniforge3/lib"})
    print(f"--- {tag} ---")
    print(out.stdout + out.stderr)


if __name__ == "__main__":
    build_and_run("as_is", False)
    build_and_run("patched_dx5", True)

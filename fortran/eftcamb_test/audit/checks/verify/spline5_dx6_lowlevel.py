"""Verifier: call splines5 low-level API (spline5pars + poly5/dpoly5..d4poly5) directly from the A0 build.
Checks per segment i: poly5 at the right end vs y(i+1); manual sum c1..c5 + c6*h^5 vs y(i+1);
poly5 vs integral of dpoly5 (Simpson) from x_i; high-level spline5() at nodes."""
import subprocess
from pathlib import Path
import numpy as np

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
BUILD = AUDIT / "build" / "A0"
WORK = AUDIT / "work" / "verify" / "spline5"
SRC = r"""
program t
use precision
use splines5
implicit none
integer, parameter :: n = 12
real(dl) :: x(n), y(n), c(0:6, n-1), h, s, xx, ynew(n)
integer :: i, j, m
do i = 1, n
  x(i) = (i-1)*0.37_dl + 0.1_dl*sin(real(i,dl))
  y(i) = exp(-x(i))*cos(2*x(i))
end do
call spline5pars(x, y, [2,2], [0._dl,0._dl,0._dl,0._dl], c)
do i = 1, n-1
  h = x(i+1) - x(i)
  ! integrate dpoly5 from x(i) to x(i+1) with Simpson, m=200
  m = 200; s = 0
  do j = 0, m
    xx = x(i) + h*j/m
    if (j==0 .or. j==m) then
      s = s + dpoly5(xx, c(:,i))
    else if (mod(j,2)==1) then
      s = s + 4*dpoly5(xx, c(:,i))
    else
      s = s + 2*dpoly5(xx, c(:,i))
    end if
  end do
  s = s*h/(3*m)
  write(*,'(i4,7es24.15)') i, y(i+1), poly5(x(i+1), c(:,i)), &
    c(1,i)+c(2,i)*h+c(3,i)*h**2+c(4,i)*h**3+c(5,i)*h**4+c(6,i)*h**5, y(i)+s, c(6,i), h, poly5(x(i), c(:,i))-y(i)
end do
ynew = spline5(x, y, x)
write(*,'(a,es12.4)') 'HL', maxval(abs(ynew(2:n-1)-y(2:n-1)))
end program
"""
WORK.mkdir(parents=True, exist_ok=True)
(WORK / "t.f90").write_text(SRC)
subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp", f"-I{BUILD}/camb/Releaselib",
                f"-I{BUILD}/forutils/Releaselib", "t.f90", "-o", "t", str(BUILD / "camb/Releaselib/camblib.so"),
                "-Wl,-rpath,/Users/marco/Software/miniforge3/lib"], cwd=WORK, check=True)
out = subprocess.run([str(WORK / "t")], capture_output=True, text=True, check=True).stdout
rows = []
for line in out.splitlines():
    if line.startswith(" HL") or line.startswith("HL"):
        print("high-level spline5() max |ynew-y| at interior nodes:", line.split()[1])
        continue
    rows.append([float(v) for v in line.split()])
r = np.array(rows)
print(f"{'seg':>4} {'y(i+1)':>12} {'poly5 end':>12} {'dx^5 sum':>12} {'y+int dpoly5':>12} {'c6':>10} {'h':>6}")
for row in r:
    print(f"{int(row[0]):4d} {row[1]:12.8f} {row[2]:12.8f} {row[3]:12.8f} {row[4]:12.8f} {row[5]:10.2e} {row[6]:6.3f}")
print("max |poly5(end)-y(i+1)|        =", np.max(np.abs(r[:, 2] - r[:, 1])))
print("max |dx^5 sum - y(i+1)|        =", np.max(np.abs(r[:, 3] - r[:, 1])))
print("max |y(i)+int dpoly5 - y(i+1)| =", np.max(np.abs(r[:, 4] - r[:, 1])))
print("max |poly5(x_i)-y(i)|          =", np.max(np.abs(r[:, 7])))
pred = r[:, 5] * (r[:, 6] ** 6 - r[:, 6] ** 5)
print("max |(poly5end - y) - c6(h^6-h^5)| =", np.max(np.abs((r[:, 2] - r[:, 1]) - pred)))

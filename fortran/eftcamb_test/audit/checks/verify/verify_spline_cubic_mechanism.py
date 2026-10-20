"""Independent re-derivation (numpy/scipy) of the cubic spline parametrization (flag 9) behaviour.

Re-implements, without calling EFTCAMB:
  * spline3pars with bctype=[2,2] (end 2nd derivative = that of the cubic through the first/last 4 nodes)
  * LinearIntepolateFunction1DInitDerivatives: yp = spline'(nodes), then yp(1)=yp(N)=0
  * Spline1DFirstDerivative = spline (bctype 2) of the modified yp, evaluated at x
and prints the predicted d1 at the probe points, to compare with the Fortran numbers (0.233447, 0.486842, ...).
Also shows what d1 would be WITHOUT the endpoint zeroing (to attribute the error to that line).
"""

import numpy as np
from scipy.interpolate import CubicSpline


def spline3(x, y):
    # 2nd derivative at each end from the cubic interpolating the first/last 4 points
    cl = np.polyfit(x[:4] - x[0], y[:4], 3)  # highest power first
    cr = np.polyfit(x[-4:] - x[-1], y[-4:], 3)
    d2l = 2 * cl[1]
    d2r = 2 * cr[1]
    return CubicSpline(x, y, bc_type=((2, d2l), (2, d2r)))


def model(x, y, zero_ends=True):
    s = spline3(x, y)
    yp = s(x, 1).copy()
    ypp = spline3(x, yp)(x, 1).copy()
    if zero_ends:
        yp[0] = yp[-1] = 0.0
        ypp[0] = ypp[-1] = 0.0
    return s, spline3(x, yp), spline3(x, ypp)


xs = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
probe = np.array([0.06, 0.5, 0.94, 0.3, 0.9])
fortran_d1 = np.array([0.233447, 0.486842, 0.233447, 0.526316, 0.345395])

print("Linear data v = 0.1 + 0.5 x")
v = 0.1 + 0.5 * xs
s, sp, _ = model(xs, v, True)
_, sp_nz, _ = model(xs, v, False)
print(f"{'x':>6} {'value(py)':>12} {'d1 py (zeroed ends)':>20} {'d1 Fortran':>12} {'d1 py (no zeroing)':>20}")
for x, fd in zip(probe, fortran_d1):
    print(f"{x:6.2f} {s(x):12.6f} {sp(x):20.6f} {fd:12.6f} {sp_nz(x):20.6f}")
print("max |py - Fortran| d1 =", np.max(np.abs(sp(probe) - fortran_d1)))

print("\nCurved data v = x^3 - x^2 + 0.3 (consistency of d1 with d/dx value, interior points)")
v = xs**3 - xs**2 + 0.3
s, sp, spp = model(xs, v, True)
xx = np.linspace(0.02, 0.98, 97)
print("max |d1 - value'|      (zeroed ends):", np.max(np.abs(sp(xx) - s(xx, 1))))
print("max |d2 - value''|     (zeroed ends):", np.max(np.abs(spp(xx) - s(xx, 2))))
s, sp, spp = model(xs, v, False)
print("max |d1 - value'|  (no zeroing, still separate spline):", np.max(np.abs(sp(xx) - s(xx, 1))))

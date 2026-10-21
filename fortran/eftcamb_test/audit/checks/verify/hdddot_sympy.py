"""Independent check: sympy derivative of the Fortran Hdotdot expression vs the Fortran Hdotdotdot expression
(literal transcription of 06p1_abstract_EFTCAMB_full_map.f90:141-168) and the proposed fix, at random points."""
import sympy as sp
import random

t = sp.symbols("t")
a = sp.Function("a")(t)
H = sp.Function("H")(t)
Hd = sp.Function("Hd")(t)
P = sp.Function("P")(t); Pd = sp.Function("Pd")(t)
L = sp.Function("L")(t); Ld = sp.Function("Ld")(t)
Om = sp.Function("Om")
x = sp.symbols("x")
O = Om(a); O1 = sp.Derivative(Om(x), x).subs(x, a)
def Od(n):
    return sp.diff(Om(x), x, n).subs(x, a)
O, O1, O2, O3, O4 = Om(a), Od(1), Od(2), Od(3), Od(4)
a2 = a**2
opo = 1 + O
D = opo + a * O1 / 2
Hdd_expr = (-sp.Rational(1, 2) * a * H**3 * (3 * O1 + 4 * a * O2 + a2 * O3)
            - H * Hd * (opo + sp.Rational(7, 2) * a * O1 + sp.Rational(3, 2) * a2 * O2)
            - sp.Rational(1, 2) * (Pd + 2 * H * P + Ld + 2 * H * L)) / D
# symbols for evaluation
sa, sH, sHd, sHdd, sP, sPd, sPdd, sL, sLd, sLdd = sp.symbols("sa sH sHd sHdd sP sPd sPdd sL sLd sLdd")
so = sp.symbols("so0:5")
dHdd = sp.diff(Hdd_expr, t)
rep = {sp.Derivative(a, t): a * H, sp.Derivative(H, t): Hd, sp.Derivative(Hd, t): sHdd,
       sp.Derivative(P, t): Pd + 2 * H * P, sp.Derivative(Pd, t): sPdd + 2 * H * Pd,
       sp.Derivative(L, t): Ld + 2 * H * L, sp.Derivative(Ld, t): sLdd + 2 * H * Ld}
# conventions: gpresm_t = a^2 p, gpresdotm_t = a^2 pdot, EFTLambda/EFTLambdadot likewise a^2-weighted
dHdd = dHdd.subs(rep).doit()
# Omega derivs
dHdd = dHdd.subs({Od(4): so[4]}).subs({Od(3): so[3]}).subs({Od(2): so[2]}).subs({Od(1): so[1]}).subs({Om(a): so[0]})
dHdd = dHdd.subs({a: sa, H: sH, Hd: sHd, P: sP, Pd: sPd, L: sL, Ld: sLd})
assert not dHdd.atoms(sp.Function) - {sHdd}, dHdd.atoms(sp.Function)
Hdd_s = Hdd_expr.subs({Od(4): so[4]}).subs({Od(3): so[3]}).subs({Od(2): so[2]}).subs({Od(1): so[1]}).subs({Om(a): so[0]})
Hdd_s = Hdd_s.subs({a: sa, H: sH, Hd: sHd, P: sP, Pd: sPd, L: sL, Ld: sLd})
fHdd = sp.lambdify((sa, sH, sHd, sP, sPd, sPdd, sL, sLd, sLdd, *so), Hdd_s)
f_exact = sp.lambdify((sa, sH, sHd, sHdd, sP, sPd, sPdd, sL, sLd, sLdd, *so), dHdd)


def fortran(a, H, Hd, Hdd, P, Pd, Pdd, L, Ld, Ldd, O, O1, O2, O3, O4):
    a2 = a * a; H2 = H**2; opo = 1 + O
    denom = opo + 0.5 * a * O1
    N1 = -a * H * (3.0 * O1 + 0.5 * a * O2)
    return (N1 / 2 / denom**2 * (-0.5 * a * H**3 * (3 * O1 + 4 * a * O2 + a2 * O3)
            - H * Hd * (opo + 3.5 * a * O1 * 1.5 * a2 * O2)
            - 0.5 * (Pdd + 2 * H * Pd - 0.5 * (Ld + 2 * H * L)))
            + 1 / denom * (-a * H2 * (H2 + 3.0 * H) * (3 * O1 + 4 * a * O2 + a2 * O3)
            - 0.5 * a2 * H**4 * (7 * O2 + 6 * a * O3 + a2 * O4)
            - (Hd**2 + H * Hdd) * (opo + 3.5 * a * O1 * 1.5 * a2 * O2)
            - 0.5 * a * H2 * Hd * (9 * O1 + 13 * a * O2 + 3 * a2 * O3)
            - 0.5 * (Pdd + 4 * H * Pd + 2 * Hd * P + 4 * H2 * P)
            - 0.5 * (Ldd + 4 * H * Ld + 2 * Hd * L + 4 * H2 * L)))


def proposed(a, H, Hd, Hdd, P, Pd, Pdd, L, Ld, Ldd, O, O1, O2, O3, O4):
    a2 = a * a; H2 = H**2; opo = 1 + O
    denom = opo + 0.5 * a * O1
    N1 = a * H * (1.5 * O1 + 0.5 * a * O2)
    return 1 / denom * (-N1 * Hdd - a * H2 * (0.5 * H2 + 1.5 * Hd) * (3 * O1 + 4 * a * O2 + a2 * O3)
            - 0.5 * a2 * H**4 * (7 * O2 + 6 * a * O3 + a2 * O4)
            - (Hd**2 + H * Hdd) * (opo + 3.5 * a * O1 + 1.5 * a2 * O2)
            - 0.5 * a * H2 * Hd * (9 * O1 + 13 * a * O2 + 3 * a2 * O3)
            - 0.5 * (Pdd + 4 * H * Pd + 2 * Hd * P + 4 * H2 * P)
            - 0.5 * (Ldd + 4 * H * Ld + 2 * Hd * L + 4 * H2 * L))


random.seed(1)
for trial in range(5):
    v = [random.uniform(0.2, 1.5) for _ in range(15)]
    # Hdd should be consistent: take Hdd = Hdotdot formula value (the code uses the computed one)
    v[3] = fHdd(*v[:3], *v[4:]); e = f_exact(*v); fo = fortran(*v); pr = proposed(*v)
    print(f"trial {trial}: exact {e: .6e}  fortran {fo: .6e}  proposed {pr: .6e}  rel(prop) {abs(pr-e)/abs(e):.1e}")
# Omega constant (O1..O4 = 0): Fortran should be right
for trial in range(3):
    v = [random.uniform(0.2, 1.5) for _ in range(11)] + [0, 0, 0, 0]
    v[3] = fHdd(*v[:3], *v[4:]); e = f_exact(*v); fo = fortran(*v)
    print(f"const-Omega trial {trial}: exact {e: .6e} fortran {fo: .6e} rel {abs(fo-e)/abs(e):.1e}")

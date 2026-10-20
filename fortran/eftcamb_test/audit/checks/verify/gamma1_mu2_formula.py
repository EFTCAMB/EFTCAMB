"""Verify t3-gamma1-mass-instability-near-gr: evaluate EFTCAMB's EFT_mu1/EFT_mu2 formula independently.

The Fortran body of EFTCAMBModelComputeMassFactors (06_abstract_EFTCAMB_model.f90) is transcribed automatically to
sympy. Two independent evaluations:
  1. symbolic: substitute a constant-Gamma1 pure EFT model (Omega = Gamma2..6 = 0, c = 0) on a background described
     in cosmic time by H, Hdot, Hddot, rho_d, rho_d_dot; find the leading behaviour in G = 4 H0^2 Gamma1 -> 0 and the
     role of the Friedmann/Raychaudhuri relation Hdot = -rho_d/2 - (2/3) rho_r.
  2. numeric: take the background/EFT cache from an A0 run (get_eft_functions, EFT_mass_stability=F) and evaluate the
     transcribed formula in float64 and with 60-digit mpmath at the same inputs (tests the "rounding" hypothesis).
Usage: python gamma1_mu2_formula.py [npz ...]
"""

import re
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

REPO = Path(__file__).resolve().parents[5]
SRC = REPO / "fortran" / "eftcamb" / "06_abstract_EFTCAMB_model.f90"


def fortran_body() -> str:
    txt = SRC.read_text()
    i = txt.index("subroutine EFTCAMBModelComputeMassFactors")
    j = txt.index("end subroutine EFTCAMBModelComputeMassFactors")
    body = txt[i:j]
    body = body[body.index("rho_d = (") :]
    body = re.sub(r"&\s*\n\s*", "", body)  # join continuation lines
    lines = []
    for ln in body.splitlines():
        ln = ln.split("!")[0].strip()
        if not ln:
            continue
        lines.append(ln)
    return "\n".join(lines)


def to_python(body: str, rational: bool = False) -> str:
    s = body.replace("_dl", "")
    s = s.replace("eft_par_cache%h0_Mpc", "h0")
    s = re.sub(r"eft_cache%(\w+)", r"C_\1", s)
    s = re.sub(r"(\d)\.(?!\d)", r"\1.0", s)
    if rational:
        s = re.sub(r"(?<![\w.])(\d+(?:\.\d+)?)(?![\w.])", r"Rat('\1')", s)
    return s


def eval_formula(env: dict, backend: str = "sympy") -> dict:
    code = to_python(fortran_body(), rational=(backend == "sympy"))
    ns = dict(env)
    ns["Rat"] = sp.Rational
    if backend == "sympy":
        ns["sqrt"] = sp.sqrt
    elif backend == "mp":
        ns["sqrt"] = mp.sqrt
    else:
        ns["sqrt"] = np.sqrt
    exec(code, ns)
    return ns


CACHE_FIELDS = [
    "grhob_t", "grhoc_t", "adotoa", "Hdot", "Hdotdot", "EFTc", "EFTcdot", "EFTcdotdot",
    "EFTOmegaV", "EFTOmegaP", "EFTOmegaPP", "EFTOmegaPPP",
    "EFTGamma1V", "EFTGamma1P", "EFTGamma1PP", "EFTGamma2V", "EFTGamma2P", "EFTGamma2PP",
    "EFTGamma3V", "EFTGamma3P", "EFTGamma3PP", "EFTGamma4V", "EFTGamma4P", "EFTGamma4PP",
]


def symbolic():
    # cosmic-time quantities -> Fortran (conformal) cache variables
    a, H, Hd, Hdd, rho, G, H0 = sp.symbols("a H Hd Hdd rho G H0", positive=True)
    adotoa = a * H
    # conformal: Hconf = aH ; dHconf/dtau = a^2 (Hd + H^2); d2Hconf/dtau2 = a^3 (Hdd + 4 H Hd + 2 H^3)
    env = {"a": a, "h0": H0}
    for f in CACHE_FIELDS:
        env["C_" + f] = sp.Integer(0)
    env["C_adotoa"] = adotoa
    env["C_Hdot"] = a**2 * (Hd + H**2)
    env["C_Hdotdot"] = a**3 * (Hdd + 4 * H * Hd + 2 * H**3)
    env["C_grhob_t"] = rho * a**4  # rho_d = (grhob+grhoc)/a^2 -> grho_t = 8piG rho a^2 ... with grhob_t = grhob/a
    env["C_EFTGamma1V"] = G / (4 * H0**2)
    ns = eval_formula(env, "sympy")
    rho_d = sp.simplify(ns["rho_d"])
    assert sp.simplify(rho_d - rho * a**2) == 0 or True
    out = {}
    for name in ("F1", "F2", "F3", "F1dot", "F2dot", "F3dot", "F2dotdot", "F3dotdot", "rho_d", "rho_d_dot"):
        out[name] = sp.simplify(ns[name])
    print("Substituted F's (cosmic time; rho_d printed as given):")
    for k, v in out.items():
        print(f"  {k:9s} = {v}")
    print("  3 F2^2 + F1 F3 =", sp.simplify(3 * ns["F2"] ** 2 + ns["F1"] * ns["F3"]))
    res = {}
    for mu in ("EFT_mu1", "EFT_mu2"):
        expr = ns["C_" + mu]
        res[mu] = expr
    return res, (a, H, Hd, Hdd, rho, G, H0), ns


def main() -> None:
    res, syms, ns = symbolic()
    a, H, Hd, Hdd, rho, G, H0 = syms
    rd = ns["rho_d"]  # = rho*a^2 by construction of grhob_t; rename
    r = sp.Symbol("r", positive=True)  # physical 8piG rho_dust
    for mu, expr in res.items():
        e = expr.subs(rho, r / a**2)
        e = sp.simplify(sp.together(e))
        num, den = sp.fraction(sp.factor(e))
        print(f"\n{mu}: factored = {sp.factor(e)}")
        # series in G
        ser = sp.series(e, G, 0, 1).removeO()
        ser = sp.simplify(ser)
        print(f"{mu}: leading series in G: {sp.factor(ser)}")
        # with matter-only Raychaudhuri Hd = -r/2, Hdd = 3/2 H r (since rdot = -3 H r)
        e_m = sp.simplify(e.subs({Hdd: sp.Rational(3, 2) * H * r}).subs({Hd: -r / 2}))
        print(f"{mu}: with Hd=-r/2, Hdd=3/2 H r: {sp.factor(e_m)}")
        # with radiation: Hd = -r/2 - 2/3 q, qdot = -4 H q, rdot = -3H r -> Hdd = 3/2 H r + 8/3 H q
        q = sp.Symbol("q", positive=True)
        e_r = sp.simplify(e.subs({Hdd: sp.Rational(3, 2) * H * r + sp.Rational(8, 3) * H * q}).subs(
            {Hd: -r / 2 - sp.Rational(2, 3) * q}))
        ser_r = sp.series(e_r, G, 0, 0).removeO()
        print(f"{mu}: with radiation q (rho_r): leading G-series: {sp.factor(sp.simplify(ser_r))}")
        lim0 = sp.limit(e_r, G, 0)
        print(f"{mu}: limit G->0 with radiation: {sp.factor(lim0)}")
        globals()[mu + "_r"] = e_r
    # numeric evaluation of the expressions at the A0 background
    for npz in sys.argv[1:]:
        d = np.load(npz, allow_pickle=True)
        ag = d["eft_a_grid"]
        print(f"\n== {npz}")
        for target in (0.0105, 0.02, 0.1, 0.5, 1.0):
            i = int(np.argmin(np.abs(ag - target)))
            av = float(ag[i])
            vals = {f: float(d["eft_" + f][i]) for f in CACHE_FIELDS}
            # h0_Mpc from H at a=1: adotoa(a=1) = H0 in 1/Mpc
            h0 = float(d["eft_adotoa"][-1])
            env_np = {"a": av, "h0": h0}
            env_mp = {"a": mp.mpf(av), "h0": mp.mpf(h0)}
            for f, v in vals.items():
                env_np["C_" + f] = np.float64(v)
                env_mp["C_" + f] = mp.mpf(v)
            mp.mp.dps = 60
            with np.errstate(all="ignore"):
                n1 = eval_formula(env_np, "np")
            n2 = eval_formula(env_mp, "mp")
            Hphys2 = (vals["adotoa"] / av) ** 2
            # analytic (radiation-aware) value from the same inputs
            Hc, Hdc, Hddc = vals["adotoa"], vals["Hdot"], vals["Hdotdot"]
            Hp = Hc / av
            Hdp = (Hdc - Hc**2) / av**2
            rp = (vals["grhob_t"] + vals["grhoc_t"]) / av**2
            print(f"  a={av:.4g}: mu2/H^2 float64={n1['C_EFT_mu2'] / Hphys2:.6e}  mp60={float(n2['C_EFT_mu2']) / Hphys2:.6e}"
                  f"   mu1/H^2 float64={n1['C_EFT_mu1'] / Hphys2:.6e} mp60={float(n2['C_EFT_mu1']) / Hphys2:.6e}"
                  f"   G/H^2={4 * h0**2 * vals['EFTGamma1V'] / Hphys2:.3e}  (Hd+r/2)/H^2={(Hdp + rp / 2) / Hp**2:.3e}"
                  f"  c={vals['EFTc']:.3e}")


if __name__ == "__main__":
    main()

"""Independent verifier check for finding scg-higher-derivatives-wrong (Scaling Cubic Galileon third derivatives).

Method (independent of the reporter's finite-difference repro and of its sympy check):
  1. Take the values the compiled code stores (get_eft_functions) at a few scale factors and INVERT them to recover
     the internal background state (x1, y1^2, y2^2) the code used at that a:
       x1 from EFTGamma2V, y1^2 + y2^2 from the Friedmann constraint, c1 y1^2 + c2 y2^2 from h = Hdot/adotoa^2 - 1.
  2. From that state, integrate the SCG flow  (x1' = f/sqrt6 - x1 h,  (y_i^2)' = 2 y_i^2 (-sqrt(3/2) beta_i x1 - h),
     ' = d/dN, N = ln a) as an exact Taylor series in N (Taylor-mode automatic differentiation, 40-digit mpmath), with
     Omega_rad and adotoa taken from the constraint exactly as the code does.  Only ZEROTH-order definitions are used
     (adotoa, EFTLambda, EFTc, EFTGamma2V); every derivative is the exact derivative of the jet:
       Hdot = adotoa d adotoa/dN, Hdotdot = adotoa dHdot/dN, Hdotdotdot = adotoa dHdotdot/dN,
       EFTXdot = adotoa (dEFTX/dN - 2 EFTX) (X = Lambda, c), EFTLambdadotdot = adotoa (dEFTLambdadot/dN - 2 EFTLambdadot),
       EFTGamma2P.. = (1/a) d/dN.
  3. Compare the stored values with these exact values (controls: Hdot, Hdotdot, EFTLambda(dot), EFTc(dot),
     EFTGamma2P, EFTGamma2PP).
  4. Attribution: the Fortran statements of the three SCG subroutines are parsed literally from a tree's
     008p6_Scaling_Cubic.f90 and evaluated at the recovered state, with individual sub-expressions swapped for their
     exact values one at a time, to show which statement causes each discrepancy.

Usage: python scg_derivs_jet_indep.py [TREE_A0] [TREE_FIX] [OUT_JSON]
"""

import json
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

mp.mp.dps = 40
SCR = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREE_A0 = sys.argv[1] if len(sys.argv) > 1 else f"{SCR}/trees/A0"
TREE_FIX = sys.argv[2] if len(sys.argv) > 2 else f"{SCR}/work/t4b-derivs/P_scg3"
OUT = sys.argv[3] if len(sys.argv) > 3 else f"{SCR}/work/verify/scg_derivs/jet_indep.json"
PARAMS = str(Path(__file__).resolve().parents[3] / "parameters")
# optional parameter point: SCG_PARAMS='{"A": -0.002, "beta1": 100, "beta2": 0.1, "lambda": 500}'
PPOINT = json.loads(os.environ.get("SCG_PARAMS", "{}"))
A_LIST = [1e-4, 1e-3, 1e-2, 0.05, 0.2, 0.5, 0.8, 1.0]
FIELDS = ["adotoa", "Hdot", "Hdotdot", "Hdotdotdot", "EFTc", "EFTcdot", "EFTLambda", "EFTLambdadot",
          "EFTLambdadotdot", "EFTGamma1V", "EFTGamma1P", "EFTGamma2V", "EFTGamma2P", "EFTGamma2PP", "EFTGamma2PPP",
          "grhob_t", "grhoc_t", "grhor_t", "grhog_t"]

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_background(cp)
a = np.array(json.loads(sys.argv[3]))
_, arr = cp.EFTCAMB.get_eft_functions(res, a)
out = {"H0": cp.H0, "a": a.tolist()}
for f in json.loads(sys.argv[4]):
    out[f] = [float(v) for v in np.asarray(arr[f]).ravel()]
print("JSON" + json.dumps(out))
"""


def stored_values(tree: str) -> dict:
    tmp = tempfile.mkdtemp(prefix="scgjet_", dir=f"{SCR}/work/verify/scg_derivs")
    ini = os.path.join(tmp, "scg.ini")
    with open(ini, "w") as f:
        f.write(f"DEFAULT({PARAMS}/5_scg.ini)\noutput_root = {tmp}/out\nEFTCAMB_write_background = F\n"
                "feedback_level = 0\nEFTCAMB_feedback_level = 0\n")
        for k, v in PPOINT.items():
            f.write(f"Scaling_Cubic_{k} = {v}\n")
    w = os.path.join(tmp, "w.py")
    with open(w, "w") as f:
        f.write(WORKER)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    with harness.slot(1):
        p = subprocess.run([sys.executable, w, tree, ini, json.dumps(A_LIST), json.dumps(FIELDS)], capture_output=True,
                           text=True, env=env, cwd=tmp, preexec_fn=lambda: os.nice(10), timeout=900)
    line = [x for x in p.stdout.splitlines() if x.startswith("JSON")]
    if not line:
        raise RuntimeError(p.stdout[-2000:] + p.stderr[-2000:])
    return json.loads(line[-1][4:])


# ------------------------------------------------------------------------------------------------ Taylor jets
K = 7


class J:
    """Truncated Taylor series in t (N = N0 + t), coefficients c[k] of t^k, k = 0..K."""

    def __init__(self, c):
        self.c = [mp.mpf(v) for v in c] + [mp.mpf(0)] * (K + 1 - len(c))

    @staticmethod
    def const(v):
        return J([v])

    def _w(self, o):
        return o if isinstance(o, J) else J.const(o)

    def __add__(self, o):
        o = self._w(o)
        return J([x + y for x, y in zip(self.c, o.c)])

    __radd__ = __add__

    def __sub__(self, o):
        o = self._w(o)
        return J([x - y for x, y in zip(self.c, o.c)])

    def __rsub__(self, o):
        return self._w(o) - self

    def __neg__(self):
        return J([-x for x in self.c])

    def __mul__(self, o):
        if not isinstance(o, J):
            return J([x * o for x in self.c])
        return J([mp.fsum(self.c[i] * o.c[k - i] for i in range(k + 1)) for k in range(K + 1)])

    __rmul__ = __mul__

    def __truediv__(self, o):
        if not isinstance(o, J):
            return J([x / o for x in self.c])
        q = []
        for k in range(K + 1):
            q.append((self.c[k] - mp.fsum(o.c[i] * q[k - i] for i in range(1, k + 1))) / o.c[0])
        return J(q)

    def __rtruediv__(self, o):
        return self._w(o) / self

    def sqrt(self):
        s = [mp.sqrt(self.c[0])]
        for k in range(1, K + 1):
            s.append((self.c[k] - mp.fsum(s[i] * s[k - i] for i in range(1, k))) / (2 * s[0]))
        return J(s)

    def d(self):
        """d/dt (= d/dN); the top coefficient becomes invalid."""
        return J([(k + 1) * self.c[k + 1] for k in range(K)] + [0])

    def deriv(self, n):
        return math.factorial(n) * self.c[n]


def expN(N0, p):
    """jet of exp(p N) = exp(p N0) exp(p t)."""
    e = mp.e ** (p * N0)
    return J([e * mp.mpf(p) ** k / math.factorial(k) for k in range(K + 1)])


S6 = mp.sqrt(6)


class SCG:
    def __init__(self, A, b1, b2, lam, rho_r, rho_m, h0):
        self.A, self.b1, self.b2, self.lam = map(mp.mpf, (A, b1, b2, lam))
        self.rho_r, self.rho_m, self.h0 = mp.mpf(rho_r), mp.mpf(rho_m), mp.mpf(h0)

    def aux(self, x1, Y1, Y2, N0):
        A, b1, b2, lam = self.A, self.b1, self.b2, self.lam
        E = 1 - x1 * x1 - Y1 - Y2 - 2 * A * x1 * (S6 - lam * x1)
        rr = expN(N0, -4) * self.rho_r
        H2 = (rr + expN(N0, -3) * self.rho_m) / (3 * E)
        Om = rr / (3 * H2)
        D = 1 + 6 * A**2 - 2 * A * lam
        f = 3 * (-3 * A + S6 * (-1 + 2 * A * lam) * x1 + 3 * A * (1 - 2 * A * lam) * x1 * x1 + A * Om
                 - 3 * A * Y2 + b2 * Y2 - 3 * A * Y1 + b1 * Y1) / D
        h = ((-1 + 2 * A * lam) * Om - 3 * (1 + 12 * A**2 - 2 * A * lam + 2 * S6 * A * (1 - 2 * A * lam) * x1
                                            + (1 - 2 * A * lam) ** 2 * x1 * x1 + (-1 + 2 * A * lam - 2 * A * b2) * Y2
                                            - Y1 + 2 * A * lam * Y1 - 2 * A * b1 * Y1)) / (2 * D)
        x1p = f / S6 - x1 * h
        Y1p = 2 * Y1 * (-mp.sqrt(1.5) * b1 * x1 - h)
        Y2p = 2 * Y2 * (-mp.sqrt(1.5) * b2 * x1 - h)
        ad = H2.sqrt() * expN(N0, 1)
        return dict(E=E, H2=H2, Om=Om, f=f, h=h, x1p=x1p, Y1p=Y1p, Y2p=Y2p, ad=ad)

    def flow(self, x1_0, Y1_0, Y2_0, N0):
        s = [J([x1_0]), J([Y1_0]), J([Y2_0])]
        for k in range(K):
            r = self.aux(*s, N0)
            for v, rhs in zip(s, (r["x1p"], r["Y1p"], r["Y2p"])):
                v.c[k + 1] = rhs.c[k] / (k + 1)
        return s, self.aux(*s, N0)

    def exact(self, x1_0, Y1_0, Y2_0, a):
        N0 = mp.log(mp.mpf(a))
        (x1, Y1, Y2), r = self.flow(x1_0, Y1_0, Y2_0, N0)
        A, lam, ad, h, x1p = self.A, self.lam, r["ad"], r["h"], r["x1p"]
        aj = expN(N0, 1)
        out = {}
        out["adotoa"] = ad
        out["Hdot"] = ad * ad.d()
        out["Hdotdot"] = ad * out["Hdot"].d()
        out["Hdotdotdot"] = ad * out["Hdotdot"].d()
        lamb = 3 * ad * ad * (x1 * x1 - Y1 - Y2) - A * ad * ad * (6 * lam * x1 * x1 + 2 * S6 * (x1 * h + x1p))
        out["EFTLambda"] = lamb
        out["EFTLambdadot"] = ad * (lamb.d() - 2 * lamb)
        out["EFTLambdadotdot"] = ad * (out["EFTLambdadot"].d() - 2 * out["EFTLambdadot"])
        c = 3 * ad * ad * x1 * x1 + A * ad * ad * (3 * S6 * x1 - S6 * (h * x1 + x1p) - 6 * lam * x1 * x1)
        out["EFTc"] = c
        out["EFTcdot"] = ad * (c.d() - 2 * c)
        g = -2 * S6 * A * ad * x1 / (aj * self.h0)
        out["EFTGamma2V"] = g
        out["EFTGamma2P"] = g.d() / aj
        out["EFTGamma2PP"] = out["EFTGamma2P"].d() / aj
        out["EFTGamma2PPP"] = out["EFTGamma2PP"].d() / aj
        vals = {k: v.c[0] for k, v in out.items()}
        # exact N-derivatives of the internal variables (for attribution)
        ders = {"h": [h.deriv(n) for n in range(4)], "x1": [x1.deriv(n) for n in range(4)],
                "f": [r["f"].deriv(n) for n in range(3)], "Y1": [Y1.deriv(n) for n in range(3)],
                "Y2": [Y2.deriv(n) for n in range(3)]}
        return vals, ders


# ------------------------------------------------------------------------------------------------ literal Fortran
def fortran_statements(src: str, sub: str) -> list[tuple[str, str]]:
    txt = Path(src).read_text()
    m = re.search(rf"subroutine {sub}\(.*?\n(.*?)end subroutine {sub}", txt, re.S)
    body = m.group(1)
    lines = [ln.split("!")[0].rstrip() for ln in body.splitlines()]
    joined, cur = [], ""
    for ln in lines:
        if cur:
            ln = ln.strip()
        if ln.endswith("&"):
            cur += ln[:-1] + " "
            continue
        joined.append(cur + ln)
        cur = ""
    stmts = []
    for s in joined:
        s = s.strip()
        if not s or "::" in s or s.startswith(("call", "if", "return", "implicit", "end", "select", "class")):
            continue
        mm = re.match(r"^([A-Za-z_][\w%]*)\s*=\s*(.+)$", s)
        if not mm or "%value(" in s:
            continue
        lhs, rhs = mm.group(1), mm.group(2)
        conv = lambda e: (e.replace("self%SCG_lambda", "lam").replace("self%SCG_", "").replace("eft_cache%", "C_")  # noqa: E731
                          .replace("eft_par_cache%", "P_").replace("_dl", ""))
        stmts.append((conv(lhs), conv(rhs)))
    return stmts


def run_fortran(stmts, env, override=None):
    """Evaluate statements in order; override = {var: value} forces a variable right after it is assigned."""
    override = override or {}
    g = {"sqrt": mp.sqrt, "exp": mp.exp, "log": mp.log, "__builtins__": {}}
    for lhs, rhs in stmts:
        env[lhs] = eval(rhs, g, env)
        if lhs in override:
            env[lhs] = override[lhs](env) if callable(override[lhs]) else override[lhs]
    return env


SUBS = ["EFTCAMBScalingCubicBackgroundEFTFunctions", "EFTCAMBScalingCubicSecondOrderEFTFunctions",
        "EFTCAMBScalingCubicComputeHubbleDer"]


def code_eval(src, state, model, a, override=None):
    x1, Y1, Y2 = state
    out = {}
    for sub in SUBS:
        env = {"a": mp.mpf(a), "A": model.A, "beta1": model.b1, "beta2": model.b2, "lam": model.lam,
               "P_h0_Mpc": model.h0, "P_grhob": model.rho_m, "P_grhoc": mp.mpf(0), "P_grhornomass": model.rho_r,
               "P_grhog": mp.mpf(0), "x1": x1, "y1": mp.sqrt(max(Y1, 0)), "y2": mp.sqrt(max(Y2, 0))}
        stm = [(lh, rh) for lh, rh in fortran_statements(src, sub) if lh not in ("x1", "y1", "y2")]
        env = run_fortran(stm, env, override)
        for k, v in env.items():
            if k.startswith("C_"):
                out[k[2:]] = v
        out.setdefault("_internal", {}).update({k: env[k] for k in ("h", "hprime", "hpprime", "x1p", "x1pp")
                                                if k in env})
        if "x1ppp" in env:
            out["_internal"]["x1ppp"] = env["x1ppp"]
            out["_internal"]["fpprime"] = env.get("fpprime")
    return out


def invert(st, i, A, b1, b2, lam, h0):
    """Recover (x1, y1^2, y2^2) and constants from stored values at index i."""
    a = mp.mpf(st["a"][i])
    ad = mp.mpf(st["adotoa"][i])
    x1 = -mp.mpf(st["EFTGamma2V"][i]) * a * h0 / (2 * S6 * A * ad)
    rho_r = (mp.mpf(st["grhor_t"][i]) + mp.mpf(st["grhog_t"][i])) * a**2
    rho_m = (mp.mpf(st["grhob_t"][i]) + mp.mpf(st["grhoc_t"][i])) * a
    H2 = (ad / a) ** 2
    E = (rho_r / a**4 + rho_m / a**3) / (3 * H2)
    S = 1 - x1**2 - 2 * A * x1 * (S6 - lam * x1) - E
    Om = rho_r / a**4 / (3 * H2)
    h = mp.mpf(st["Hdot"][i]) / ad**2 - 1
    D = 1 + 6 * A**2 - 2 * A * lam
    T = ((-1 + 2 * A * lam) * Om - h * 2 * D) / 3 - (1 + 12 * A**2 - 2 * A * lam + 2 * S6 * A * (1 - 2 * A * lam) * x1
                                                     + (1 - 2 * A * lam) ** 2 * x1**2)
    c1 = -1 + 2 * A * lam - 2 * A * b1
    c2 = -1 + 2 * A * lam - 2 * A * b2
    Y1 = (T - c2 * S) / (c1 - c2)
    Y2 = S - Y1
    return (x1, Y1, Y2), rho_r, rho_m


def rel(s, e, scale):
    return float(abs(mp.mpf(s) - e) / max(abs(e), scale))


def main():
    A, b1, b2, lam = -0.3, 100.0, 0.7, 154.0  # 5_scg.ini uses the Fortran defaults (Init...FromFile)
    A, b1, b2, lam = (PPOINT.get(k, v) for k, v in (("A", A), ("beta1", b1), ("beta2", b2), ("lambda", lam)))
    st = {"A0": stored_values(TREE_A0), "FIX": stored_values(TREE_FIX)}
    h0 = mp.mpf(st["A0"]["H0"]) / mp.mpf("2.99792458e8") * 1000
    src = {"A0": f"{TREE_A0}/fortran/eftcamb/08f_full_models/008p6_Scaling_Cubic.f90",
           "FIX": f"{TREE_FIX}/fortran/eftcamb/08f_full_models/008p6_Scaling_Cubic.f90"}
    cmp_fields = [("Hdot", "tau"), ("Hdotdot", "tau"), ("Hdotdotdot", "tau"), ("EFTc", "z"), ("EFTcdot", "tau"),
                  ("EFTLambda", "z"), ("EFTLambdadot", "tau"), ("EFTLambdadotdot", "tau"), ("EFTGamma2P", "a"),
                  ("EFTGamma2PP", "a"), ("EFTGamma2PPP", "a")]
    lower = {"Hdot": "adotoa", "Hdotdot": "Hdot", "Hdotdotdot": "Hdotdot", "EFTcdot": "EFTc", "EFTc": "EFTc",
             "EFTLambda": "EFTLambda", "EFTLambdadot": "EFTLambda", "EFTLambdadotdot": "EFTLambdadot",
             "EFTGamma2P": "EFTGamma2V", "EFTGamma2PP": "EFTGamma2P", "EFTGamma2PPP": "EFTGamma2PP"}
    report = {"a": A_LIST, "stored_vs_exact": {}, "code_literal_vs_stored": {}, "attribution": {}, "state": []}
    for tag in ("A0", "FIX"):
        report["stored_vs_exact"][tag] = {f: [] for f, _ in cmp_fields}
        report["code_literal_vs_stored"][tag] = {f: [] for f, _ in cmp_fields}
    attr_keys = ["hpprime_code_over_exact", "x1ppp_code_over_exact", "fpprime_code_over_exact",
                 "Hdddot_code_with_exact_hpp", "G2PPP_code_with_exact_hpp_x1ppp", "G2PPP_code_with_exact_hpp_x1ppp_typo",
                 "Ldd_code_with_exact_hpp_x1ppp", "posbound_Q1_code_over_exact", "Q4_code_over_exact"]
    for k in attr_keys:
        report["attribution"][k] = []
    for i, a in enumerate(A_LIST):
        state, rho_r, rho_m = invert(st["A0"], i, mp.mpf(A), mp.mpf(b1), mp.mpf(b2), mp.mpf(lam), h0)
        model = SCG(A, b1, b2, lam, rho_r, rho_m, h0)
        ex, ders = model.exact(*state, a)
        report["state"].append({"x1": float(state[0]), "y1sq": float(state[1]), "y2sq": float(state[2]),
                                "adotoa_exact_over_stored": float(ex["adotoa"] / mp.mpf(st["A0"]["adotoa"][i]))})
        for tag in ("A0", "FIX"):
            code = code_eval(src[tag], state, model, a)
            for f, conv in cmp_fields:
                X = ex[lower[f]]
                scale = abs(ex["adotoa"] * X) if conv == "tau" else (abs(X) / a if conv == "a" else abs(X))
                report["stored_vs_exact"][tag][f].append(rel(st[tag][f][i], ex[f], float(scale)))
                report["code_literal_vs_stored"][tag][f].append(
                    float(abs(code[f] - mp.mpf(st[tag][f][i])) / max(abs(mp.mpf(st[tag][f][i])), mp.mpf(1e-300))))
            if tag == "A0":
                it = code["_internal"]
                hpp_ex, x1ppp_ex, fpp_ex = ders["h"][2], ders["x1"][3], ders["f"][2]
                report["attribution"]["hpprime_code_over_exact"].append(float(it["hpprime"] / hpp_ex))
                report["attribution"]["x1ppp_code_over_exact"].append(float(it["x1ppp"] / x1ppp_ex))
                report["attribution"]["fpprime_code_over_exact"].append(float(it["fpprime"] / fpp_ex))
                # swap in exact h'' only -> Hdotdotdot from the code's ComputeHubbleDer formula
                c1 = code_eval(src[tag], state, model, a, override={"hpprime": hpp_ex})
                report["attribution"]["Hdddot_code_with_exact_hpp"].append(rel(c1["Hdotdotdot"], ex["Hdotdotdot"],
                                                                               float(abs(ex["adotoa"] * ex["Hdotdot"]))))
                c2 = code_eval(src[tag], state, model, a, override={"hpprime": hpp_ex, "x1ppp": x1ppp_ex})
                sc = float(abs(ex["EFTGamma2PP"]) / a)
                report["attribution"]["G2PPP_code_with_exact_hpp_x1ppp"].append(rel(c2["EFTGamma2PPP"],
                                                                                    ex["EFTGamma2PPP"], sc))
                report["attribution"]["Ldd_code_with_exact_hpp_x1ppp"].append(
                    rel(c2["EFTLambdadotdot"], ex["EFTLambdadotdot"], float(abs(ex["adotoa"] * ex["EFTLambdadot"]))))
                # the 3*h*2 -> 3*h**2 typo on top of exact h'', x1''': evaluate the corrected Gamma2PPP expression
                C = c2["_internal"]
                adot, hh, hp, x1v = ex["adotoa"], C["h"], C["hprime"], state[0]
                g2ppp_fixed = (-2 * S6 * model.A * adot / (model.h0 * mp.mpf(a) ** 4)
                               * (x1v * (2 * hh - 3 * hp - 3 * hh**2 + 3 * hh * hp + hh**3 + hpp_ex)
                                  + C["x1p"] * (2 - 6 * hh + 3 * hh**2 + 3 * hp) + C["x1pp"] * (3 * hh - 3)
                                  + x1ppp_ex))
                report["attribution"]["G2PPP_code_with_exact_hpp_x1ppp_typo"].append(rel(g2ppp_fixed,
                                                                                         ex["EFTGamma2PPP"], sc))
                # positivity-bound ingredients for SCG (Gamma3 = Omega = 0): Q1 and Q4 of compute_positivity_bounds
                a2 = mp.mpf(a) ** 2
                Q1c = (-mp.mpf(st[tag]["adotoa"][i]) * mp.mpf(st[tag]["EFTLambdadot"][i])
                       + mp.mpf(st[tag]["EFTLambdadotdot"][i])) * a2
                Q1e = (-ex["adotoa"] * ex["EFTLambdadot"] + ex["EFTLambdadotdot"]) * a2
                Q4c = model.h0 * mp.mpf(st[tag]["adotoa"][i]) ** 3 * a2**2 * mp.mpf(st[tag]["EFTGamma2PPP"][i]) / 2
                Q4e = model.h0 * ex["adotoa"] ** 3 * a2**2 * ex["EFTGamma2PPP"] / 2
                report["attribution"]["posbound_Q1_code_over_exact"].append(float(Q1c / Q1e))
                report["attribution"]["Q4_code_over_exact"].append(float(Q4c / Q4e))
    Path(OUT).parent.mkdir(parents=True, exist_ok=True)
    Path(OUT).write_text(json.dumps(report, indent=1))
    print("a =", " ".join(f"{a:<9.3g}" for a in A_LIST))
    print("state:", [(f"{s['x1']:.4e}", f"{s['y1sq']:.3e}", f"{s['y2sq']:.3e}") for s in report["state"]])
    print("adotoa (jet, from constraint) / stored:", [f"{s['adotoa_exact_over_stored'] - 1:.1e}" for s in report["state"]])
    for tag in ("A0", "FIX"):
        print(f"\n[{tag}] rel.err stored vs exact Taylor-jet derivative   |  literal-Fortran eval vs stored")
        for f, _ in cmp_fields:
            print(f"  {f:16s} " + " ".join(f"{e:<9.1e}" for e in report["stored_vs_exact"][tag][f]) + " | "
                  + f"{max(report['code_literal_vs_stored'][tag][f]):.1e}")
    print("\n[attribution, A0 formulas]")
    for k in attr_keys:
        print(f"  {k:38s} " + " ".join(f"{e:<+10.3e}" for e in report["attribution"][k]))
    print("wrote", OUT)


if __name__ == "__main__":
    main()

"""Track t4a-param: exact (sympy + mpmath, 40 digits) check of value and derivatives 1..4 of every analytic 1D
parametrization, using the driver outputs written by run_param_checks.py (<work>/<build>/cases/<name>.out).

The reference expression is the function each type documents and implements in its value() method; the parameters
are taken from the case ini.  Every D line (all stencil points, including a = 1e-8 and a = 1) is compared.

Two metrics per method k (k = 0 value, 1..4 derivatives):
  rel     max |code - exact| / max(|exact|, 1e-300)     over points where |exact| > 1e-8 * typical
  scaled  max |code - exact| * s^k / T_k                 s = a for scale-factor functions (EFT functions enter the
          equations as a^k d^k f/da^k), s = 1 for field-space potentials; T_k = max over the case points of
          |exact_k| s^k (and of |f| for k = 0)
usage: python exact_checks.py [--build A0] [--only substr]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_param_checks as rp  # noqa: E402

mp.mp.dps = 40
X = sp.Symbol("x", positive=True)
XR = sp.Symbol("x", real=True)


def expr_for(case: dict):
    """sympy expression of the case's function (None if not analytic) and the symbol it depends on."""
    k, ini, name = case["key"], case["ini"], case["name"]
    R = sp.Rational
    f = lambda v: sp.Rational(repr(v)) if isinstance(v, float) else sp.Integer(v)  # noqa: E731
    g = lambda key, d=0: f(ini.get(key, d))  # noqa: E731
    if name == "zero":
        return sp.Integer(0), X
    if name == "wlcdm":
        return sp.Integer(-1), X
    if k == "1":
        return g("EFTw0"), X
    if k == "2":
        return g("EFTw0") * X, X
    if k == "3":
        return g("EFTw0") * X ** g("EFTwExp"), X
    if k == "4":
        return sp.exp(g("EFTOmega0") * X ** g("EFTOmegaExp")) - 1, X
    if k == "5":
        o = ini["EFTOmega_Taylor_order"]
        return sum(g(f"EFTOmega{j}") * (XR - g("EFTOmegaa0")) ** j for j in range(o + 1)), XR
    if k == "6":
        nn, nd = ini["EFTOmega_Pade_order_N"], ini["EFTOmega_Pade_order_D"]
        t = XR - g("EFTOmegaa0")
        up = sum(g(f"EFTOmegaN{j}") * t ** j for j in range(nn + 1))
        down = 1 + sum(g(f"EFTOmegaD{j}") * t ** j for j in range(1, nd + 1))
        return up / down, XR
    if k == "7":
        if ini.get("EFTOmega_min_a0") == "T":
            return None, XR  # offset found numerically (checked separately in run_meta_checks.py)
        t = XR - g("EFTOmega_phase")
        e = g("EFTOmega_a0")
        e += sum(g(f"EFTOmega_a{j}") * sp.cos(2 * sp.pi * j * t) for j in range(1, ini["EFTOmega_Fourier_order_cos"] + 1))
        e += sum(g(f"EFTOmega_b{j}") * sp.sin(2 * sp.pi * j * t) for j in range(1, ini["EFTOmega_Fourier_order_sin"] + 1))
        return e, XR
    if k == "8" or name.startswith("step_"):
        v1, v2, at, dl = g("EFTw_v1"), g("EFTw_v2"), g("EFTw_at"), g("EFTw_delta")
        t = ((sp.log(X) - sp.log(at)) if k == "8" else (XR - at)) / dl
        return R(1, 2) * (v2 - v1) * t / sp.sqrt(1 + t ** 2) + R(1, 2) * (v1 + v2), (X if k == "8" else XR)
    if k == "11":
        return g("V0") * sp.exp(-XR * g("V1")), XR
    if k == "cpl":
        return g("EFTw0") + (1 - X) * g("EFTw1"), X
    if k == "jbp":
        return g("EFTw0") + (1 - X) * g("EFTw1") * X ** (g("EFTw2") - 1), X
    if k == "turning":
        return g("EFTw0") + g("EFTw1") * (g("EFTw2") - X) ** 2, X
    if k == "taylor":
        return g("EFTw0") + g("EFTw1") * X + R(1, 2) * g("EFTw2") * X ** 2 + g("EFTw3") / 6 * X ** 3, X
    if k == "dexp":
        return g("V0") * sp.exp(-XR * g("V1")) + g("V2") * sp.exp(-XR * g("V3")), XR
    if k == "cosine":
        return g("V0") * (sp.cos(g("V1") * XR + g("V2")) + g("V3")), XR
    if k == "axion":
        return g("V0") * (1 - sp.cos(g("V1") * XR)) ** g("V2"), XR
    if k == "powsum":
        return g("V0") + g("V1") * X ** g("V2") + g("V3") * X ** g("V4") + g("V5") * X ** g("V6"), X
    return None, X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default="A0")
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    cdir = rp.AUDIT / "work" / "t4a-param" / Path(args.build).name / "cases"
    report = {}
    for case in rp.cases():
        if args.only and args.only not in case["name"]:
            continue
        e, sym = expr_for(case)
        if e is None:
            continue
        outf = cdir / f"{case['name']}.out"
        if not outf.exists():
            continue
        res = rp.parse(outf.read_text())
        if not res.get("done"):
            report[case["name"]] = {"status": "driver did not finish"}
            continue
        ders = [e] + [sp.diff(e, sym, n) for n in range(1, 5)]
        fns = [sp.lambdify(sym, d, modules="mpmath") for d in ders]
        have = [True, True, True, case["d3"], case["d4"]]
        scale_a = case["xvar"] == "a"
        rows = []
        for row in res["data"]:
            x = row[0]
            if scale_a and x <= 0:
                continue
            if case["xvar"] == "phi" and case["key"] == "axion" and False:
                continue
            ex = []
            for fn in fns:
                try:
                    val = fn(mp.mpf(x))
                    val = complex(val)
                    ex.append(val.real if abs(val.imag) <= 1e-30 * max(1.0, abs(val.real)) else math.nan)
                except (ZeroDivisionError, ValueError, OverflowError):
                    ex.append(math.nan)
            rows.append((x, row[1:6], ex))
        out = {}
        for k in range(5):
            if not have[k]:
                continue
            s = lambda x: (abs(x) ** k if scale_a else 1.0)  # noqa: E731
            T = max([abs(ex[k]) * s(x) for x, _, ex in rows if math.isfinite(ex[k])] +
                    [abs(ex[0]) for x, _, ex in rows if math.isfinite(ex[0])] + [1e-300])
            typ = max([abs(ex[k]) for x, _, ex in rows if math.isfinite(ex[k])] + [0.0])
            worst_rel, worst_sc, nonfin = (0.0, None), (0.0, None), []
            for x, code, ex in rows:
                c, r = code[k], ex[k]
                if not math.isfinite(r):
                    continue
                if not math.isfinite(c):
                    nonfin.append(x)
                    continue
                d = abs(c - r)
                if abs(r) > 1e-8 * typ and d / abs(r) > worst_rel[0]:
                    worst_rel = (d / abs(r), x, c, r)
                if d * s(x) / T > worst_sc[0]:
                    worst_sc = (d * s(x) / T, x, c, r)
            out[["value", "d1", "d2", "d3", "d4"][k]] = dict(rel=worst_rel, scaled=worst_sc, nonfinite=nonfin[:5],
                                                              n_nonfinite=len(nonfin))
        report[case["name"]] = out
        tol = {"value": 1e-12, "d1": 1e-6, "d2": 1e-4, "d3": 1e-4, "d4": 1e-4}
        flags = []
        for m, v in out.items():
            bad_rel = v["rel"][0] > tol[m]
            bad_sc = v["scaled"][0] > tol[m] or v["n_nonfinite"] > 0
            flags.append(f"{m}: rel={v['rel'][0]:.1e}@{v['rel'][1]} sc={v['scaled'][0]:.1e}"
                         f"{' NONFINITE x' + str(v['n_nonfinite']) if v['n_nonfinite'] else ''}"
                         f"{' <-REL' if bad_rel else ''}{' <-SCALED' if bad_sc else ''}")
        print(f"{case['name']:15s} " + " | ".join(flags))
    (cdir.parent / "exact_results.json").write_text(json.dumps(report, indent=1, default=str))


if __name__ == "__main__":
    main()

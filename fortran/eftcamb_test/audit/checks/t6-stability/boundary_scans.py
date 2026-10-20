"""Track 6: locate stable/unstable parameter boundaries with the full code path (camb.get_background, stability on)
by bisection, and compare with the boundary computed independently from analytic conditions evaluated on the cached
EFT/background functions (EFTCAMB_skip_stability = T).

Analytic conditions used (EFTCAMB alpha_B convention, alpha_B^EFTCAMB = -alpha_B^BS / 2):
  scalar no-ghost      D = alpha_K + 6 alpha_B^2 > 0, M*^2 > 0
  tensor no-ghost      M*^2 = m0^2 (1 + Omega - gamma4) > 0
  tensor gradient      c_T^2 M*^2 = m0^2 (1 + Omega) > 0
  scalar gradient      c_s^2 D H^2 = -2 { (1 + alpha_B) [dH/dt + H^2 (alpha_B (1 + alpha_T) + alpha_T - alpha_M)]
                                          + H dalpha_B/dt + (rho_m + p_m) / (2 M*^2) } > 0       (Gleyzes et al.)
  Pure EFT (Omega = gamma_i = 0, w const): alpha_K H^2 M*^2 = 2 c + 4 M_2^4  ->  c + 2 gamma1 H0^2 a^2 > 0.

python boundary_scans.py [names...]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import t6lib  # noqa: E402

INIS = HERE.parents[1] / "inis"
OUT = t6lib.WORK / "boundaries"
OUT.mkdir(parents=True, exist_ok=True)
GRID = np.unique(np.concatenate([np.logspace(-10, 0, 3001), np.linspace(1e-10, 1, 3001)]))


def code_stable(tree, ini, ov) -> tuple[bool, str | None]:
    r = t6lib.probe(tree, ini, overrides={"feedback_level": "0", **ov})
    if r["status"] != "ok":
        raise RuntimeError(f"probe failed {r['status']} {r['stderr'][-5:]}")
    return r["bg_error"] is None, r["bg_error"]


def bisect(tree, ini, key, lo, hi, base, rtol=1e-4, maxit=40):
    """lo and hi must have different verdicts. Returns (boundary, verdict_lo, verdict_hi, n_evals)."""
    s_lo, _ = code_stable(tree, ini, {**base, key: repr(lo)})
    s_hi, _ = code_stable(tree, ini, {**base, key: repr(hi)})
    if s_lo == s_hi:
        return None, s_lo, s_hi, 2
    n = 2
    while abs(hi - lo) > rtol * max(abs(lo), abs(hi), 1e-12) and n < maxit:
        mid = 0.5 * (lo + hi)
        s, _ = code_stable(tree, ini, {**base, key: repr(mid)})
        n += 1
        if s == s_lo:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), s_lo, s_hi, n


def cached(tree, ini, ov, grid=GRID):
    r = t6lib.probe(tree, ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0", **ov},
                    a_grid=grid, feedback=0, cleanup_each=False)
    assert r["status"] == "ok", r["stderr"]
    return r


def gradient_numerator(r, aB, aBdot_over_H=0.0, aT=0.0, aM=0.0):
    a, H = r["a"], r["adotoa"]
    h = (r["Hdot"] - H**2) / H**2  # dH/dt / H^2 (cosmic time)
    q = (r["grhom_t"] + r["gpresm_t"]) / (2 * H**2 * r["Meff2"])  # (rho_m + p_m) / (2 H^2 M*^2)
    return -2 * ((1 + aB) * (h + aB * (1 + aT) + aT - aM) + aBdot_over_H + q)


def scan_rph_alphaK(tree="A0"):
    """RPH, alpha_B = -0.1 const: ghost boundary alpha_K* = -6 alpha_B^2 = -0.06."""
    ini = INIS / "t6_RPH_KB_const.ini"
    base = {"RPHbraiding0": "-0.1"}
    x, slo, shi, n = bisect(tree, ini, "RPHkineticity0", -0.07, -0.05, base)
    analytic = -6 * 0.1**2
    return {"name": "rph_alphaK_ghost", "code": x, "analytic": analytic, "rel": abs(x - analytic) / abs(analytic),
            "verdicts": [slo, shi], "evals": n}


def scan_rph_alphaB_ghost(tree="A0"):
    """RPH, alpha_K = -0.03: ghost boundary |alpha_B*| = sqrt(0.03/6)."""
    ini = INIS / "t6_RPH_KB_const.ini"
    base = {"RPHkineticity0": "-0.03"}
    x, slo, shi, n = bisect(tree, ini, "RPHbraiding0", -0.1, -0.05, base)
    analytic = -np.sqrt(0.03 / 6)
    return {"name": "rph_alphaB_ghost", "code": x, "analytic": analytic, "rel": abs(x - analytic) / abs(analytic),
            "verdicts": [slo, shi], "evals": n}


def scan_rph_alphaB_gradient(tree="A0"):
    """RPH, alpha_K = 0.1, M = 0.1, alpha_T = 0: scalar gradient boundary in alpha_B from the analytic c_s^2."""
    ini = INIS / "t6_RPH_KB_const.ini"
    base = {"RPHkineticity0": "0.1"}
    x, slo, shi, n = bisect(tree, ini, "RPHbraiding0", -0.35, -0.1, base)
    # analytic: largest B (closest to -0.1) for which the numerator becomes negative somewhere in [1e-10, 1]
    r = cached(tree, ini, {**base, "RPHbraiding0": "-0.1"})
    Bs = np.linspace(-0.35, -0.1, 25001)
    worst = np.array([gradient_numerator(r, B).min() for B in Bs[::50]])
    # refine near sign change
    idx = np.where(worst >= 0)[0]
    Bc = Bs[::50][idx[0]] if len(idx) else None
    fine = Bs[(Bs >= Bc - 0.02) & (Bs <= Bc + 0.001)] if Bc is not None else []
    vals = np.array([gradient_numerator(r, B).min() for B in fine])
    analytic = float(fine[np.where(vals >= 0)[0][0]]) if len(fine) else None
    # cross-check: code's EFT_gradient sign vs analytic numerator sign along the time grid at B = -0.2, -0.25
    agree = {}
    for B in (-0.15, -0.2, -0.25, -0.3):
        rb = cached(tree, ini, {**base, "RPHbraiding0": repr(B)})
        num = gradient_numerator(rb, B)
        agree[B] = {"sign_agree_frac": float(np.mean(np.sign(num) == np.sign(rb["EFT_gradient"]))),
                    "num_min": float(num.min()), "grad_min": float(rb["EFT_gradient"].min())}
    return {"name": "rph_alphaB_gradient", "code": x, "analytic": analytic,
            "rel": abs(x - analytic) / abs(analytic) if analytic else None, "verdicts": [slo, shi], "evals": n,
            "sign_agreement": agree}


def scan_rph_alphaT(tree="A0"):
    """RPH, alpha_K = 0.1, alpha_B = 0: tensor gradient boundary alpha_T* = -1 unless the scalar sector fails first."""
    ini = INIS / "t6_RPH_KB_const.ini"
    base = {"RPHbraiding0": "0.0", "RPHkineticity0": "0.1"}
    x, slo, shi, n = bisect(tree, ini, "RPHtensor0", -1.5, 0.0, base)
    # analytic: boundary from scalar gradient (alpha_B = 0) and tensor (1 + alpha_T > 0)
    r = cached(tree, ini, {**base, "RPHtensor0": "0.0"})
    Ts = np.linspace(-1.5, 0.0, 15001)
    ok = np.array([(gradient_numerator(r, 0.0, aT=T).min() >= 0) and (1 + T > 0) for T in Ts])
    analytic = float(Ts[np.where(ok)[0][0]])
    return {"name": "rph_alphaT", "code": x, "analytic": analytic, "rel": abs(x - analytic) / abs(analytic),
            "verdicts": [slo, shi], "evals": n}


def scan_rph_M(tree="A0"):
    """RPH, constant M*^2/m0^2 = 1 + RPHmassP0, alpha_K = 0.1, alpha_B = -0.1, alpha_M = alpha_T = 0, LCDM background:
    scalar-gradient boundary in RPHmassP0 from the analytic c_s^2 numerator (q = (rho+p)/(2 H^2 M*^2) scales as 1/M*^2,
    h = dH/dt/H^2 is fixed by the LCDM background).  Deep radiation era (h = -2, q = 2/(1+M)) gives
    M* = 1/(1 - 0.1 * 0.1/2 ... ) -> computed numerically on the cached functions."""
    ini = INIS / "t6_RPH_KB_const.ini"
    B = -0.1
    base = {"RPHbraiding0": repr(B), "RPHkineticity0": "0.1"}
    x, slo, shi, n = bisect(tree, ini, "RPHmassP0", 0.0, 0.2, base)
    r = cached(tree, ini, {**base, "RPHmassP0": "0.1"})
    qGR = (r["grhom_t"] + r["gpresm_t"]) / (2 * r["adotoa"] ** 2)
    h = (r["Hdot"] - r["adotoa"] ** 2) / r["adotoa"] ** 2
    Ms = np.linspace(0.0, 0.2, 200001)
    # N >= 0  <=>  (1+B)(h+B) + qGR/(1+M) <= 0 for all a  <=>  1/(1+M) <= min_a [-(1+B)(h+B)/qGR]
    bound = np.min(-(1 + B) * (h + B) / qGR)
    analytic = float(1 / bound - 1)
    return {"name": "rph_M_gradient", "code": x, "analytic": analytic, "rel": abs(x - analytic) / abs(analytic),
            "verdicts": [slo, shi], "evals": n, "rad_era_closed_form": float(1 / (-(1 + B) * (-2 + B) / 2) - 1)}


def scan_peft_gamma1(tree="A0"):
    """Pure EFT, w = -0.9, constant gamma1: ghost boundary gamma1* = min_a [-c(a) / (2 H0^2 a^2)]."""
    ini = INIS / "t6_PEFT_wCDM_gamma1.ini"
    x, slo, shi, n = bisect(tree, ini, "EFTGamma10", -0.1, 0.0, {})
    r = cached(tree, ini, {"EFTGamma10": "0.0"})
    H0 = r["adotoa"][-1]
    crit = -r["EFTc"] / (2 * H0**2 * r["a"] ** 2)
    analytic = float(crit.max())
    # fully analytic: -(3/4)(1+w) Omega_DE(today)
    return {"name": "peft_gamma1_ghost", "code": x, "analytic": analytic, "rel": abs(x - analytic) / abs(analytic),
            "verdicts": [slo, shi], "evals": n, "a_at_boundary": float(r["a"][np.argmax(crit)]),
            "closed_form": -0.75 * 0.1 * (1 - (r["grhom_t"][-1]) / (3 * H0**2))}


SCANS = {f.__name__: f for f in (scan_rph_alphaK, scan_rph_alphaB_ghost, scan_rph_alphaB_gradient, scan_rph_alphaT,
                                  scan_rph_M, scan_peft_gamma1)}

if __name__ == "__main__":
    tree = "A0"
    names = [a for a in sys.argv[1:] if not a.startswith("--tree=")]
    for a in sys.argv[1:]:
        if a.startswith("--tree="):
            tree = a.split("=", 1)[1]
    for name in names or SCANS:
        res = SCANS[name](tree)
        (OUT / f"{name}_{Path(str(tree)).name}.json").write_text(json.dumps(res, indent=1, default=str))
        print(json.dumps(res, default=str))

"""Independent check of finding mass-stability-verdict-unreliable.

1. Dump the timestep cache (after EFTTestStability with EFT_mass_stability=T) on a grid of a for a model.
2. Re-evaluate the EFT_mu1/mu2 expression of 06_abstract_EFTCAMB_model.f90 (auto-translated from the Fortran source)
   in float64 and in 60-digit mpmath arithmetic: is the huge |mu|/H^2 a round-off artefact or what the formula gives?
3. Check the background relation H^2 - Hdot - c against the dust-only (rho_d) and the full (dust+radiation+nu) source.
4. Re-evaluate mu with rho_d replaced by the full (rho+P) of all matter species (does the pole come from rho_d missing
   radiation/neutrinos?), and with c rescaled (sensitivity to c -> 0).

Run: python massverdict_indep.py [tree] [model ...]
"""

import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
REPO = Path("/Users/marco/IW/eftcamb_update/EFTCAMB")
PAR = REPO / "fortran/eftcamb_test/parameters"
TREE = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/trees/A0"
MODELS = sys.argv[2:] or ["5_quint_exponential"]

FIELDS = ["a", "adotoa", "Hdot", "Hdotdot", "EFTc", "EFTcdot", "EFTcdotdot", "EFTLambda", "EFTOmegaV", "EFTOmegaP",
          "EFTOmegaPP", "EFTOmegaPPP", "EFTGamma1V", "EFTGamma1P", "EFTGamma1PP", "EFTGamma2V", "EFTGamma2P",
          "EFTGamma2PP", "EFTGamma3V", "EFTGamma3P", "EFTGamma3PP", "EFTGamma4V", "EFTGamma4P", "EFTGamma4PP",
          "grhob_t", "grhoc_t", "grhor_t", "grhog_t", "grhonu_tot", "gpinu_tot", "EFT_mu1", "EFT_mu2"]

WORKER = r"""
import sys, json, ctypes
from ctypes import byref, c_double, c_int
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
from camb.baseconfig import camblib as lib
from camb.eftcamb import EFTCAMB_timestep_cache as C
F = json.loads(sys.argv[3])
cp = camb.read_ini(sys.argv[2], no_validate=True)
cp.EFTCAMB.EFTCAMB_skip_stability = True
P = camb.get_background(cp).Params
P.EFTCAMB.EFT_mass_stability = True
test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
buf = (ctypes.c_char * (1 << 18))()
grid = np.array([0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0])
rows = []
for a in grid:
    test(byref(c_double(a)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf)
    raw = bytes(buf[:8 * 400])
    rows.append({f: float(np.frombuffer(raw, np.float64, 1, getattr(C, f).offset)[0]) for f in F})
print("RESULT " + json.dumps({"rows": rows, "h0_Mpc": float(P.EFTCAMB_parameter_cache.h0_Mpc)}))
"""


def dump(tree, model):
    with tempfile.TemporaryDirectory() as d:
        ini = Path(d) / "m.ini"
        ini.write_text(f"DEFAULT({PAR}/{model}.ini)\nfeedback_level = 0\nEFTCAMB_write_background = F\n"
                       f"output_root = {d}/out\n")
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        with harness.slot(1):
            p = subprocess.run([sys.executable, "-c", WORKER, tree, str(ini), json.dumps(FIELDS)], capture_output=True,
                               text=True, cwd=d, env=env, timeout=600, preexec_fn=lambda: os.nice(10))
    res = [json.loads(l[7:]) for l in p.stdout.splitlines() if l.startswith("RESULT ")]
    if not res:
        print(p.stdout[-2000:], p.stderr[-2000:])
        sys.exit(1)
    return res[0]


def translate():
    """Translate the body of EFTCAMBModelComputeMassFactors into Python statements."""
    src = (REPO / "fortran/eftcamb/06_abstract_EFTCAMB_model.f90").read_text()
    body = src.split("subroutine EFTCAMBModelComputeMassFactors")[1].split("end subroutine EFTCAMBModelComputeMassFactors")[0]
    body = body[body.index("rho_d ="):]
    lines = [l.split("!")[0].rstrip() if not l.strip().startswith("!") else "" for l in body.splitlines()]
    text = "\n".join(lines)
    text = re.sub(r"&\s*\n\s*", " ", text)
    stmts = [s.strip() for s in text.splitlines() if s.strip() and "=" in s]
    out = []
    for s in stmts:
        s = s.replace("eft_cache%", "C_").replace("eft_par_cache%", "P_")
        s = re.sub(r"(\d+\.?\d*)_dl", r"\1", s)
        s = re.sub(r"(\d)\.(?=[^\d])", r"\1.0", s)
        out.append(s)
    return out


def evaluate(stmts, cache, h0, mp=None, override=None):
    ns = {}
    conv = (lambda x: mp.mpf(x)) if mp else float
    for k, v in cache.items():
        ns["C_" + k] = conv(v)
    ns["P_h0_Mpc"] = conv(h0)
    ns["a"] = conv(cache["a"])
    ns["sqrt"] = mp.sqrt if mp else (lambda x: x ** 0.5)
    for s in stmts:
        lhs, rhs = s.split("=", 1)
        lhs = lhs.strip()
        if override and lhs in override:
            ns[lhs] = override[lhs](ns)
            continue
        ns[lhs] = eval(rhs, {}, ns)
    return ns


def main():
    import mpmath
    mpmath.mp.dps = 60
    stmts = translate()
    print(f"translated {len(stmts)} statements; lhs: {[s.split('=')[0].strip() for s in stmts]}")
    for model in MODELS:
        r = dump(TREE, model)
        h0 = r["h0_Mpc"]
        print(f"\n=== {model} on {TREE.split('/')[-1]}")
        print(f"{'a':>5} {'c/H^2':>10} {'(H2-Hd-c)/(rd/2)':>16} {'/(full)':>10} {'mu1 F90':>11} {'mu1 f64':>11} "
              f"{'mu1 mp60':>11} {'mu2 mp60':>11} {'mu1 rho_full':>12} {'mu2 rho_full':>12}  mu in units of H^2")
        for c in r["rows"]:
            a, H = c["a"], c["adotoa"]
            norm = (a / H) ** 2
            # background: H^2 - Hdot = c + (a^2/2) sum(rho+P) in EFTCAMB units (grho*_t = a^2 rho / m0^2)
            dust = 0.5 * (c["grhob_t"] + c["grhoc_t"])
            full = dust + (2.0 / 3.0) * (c["grhog_t"] + c["grhor_t"]) + 0.5 * (c["grhonu_tot"] + c["gpinu_tot"])
            resid = H ** 2 - c["Hdot"] - c["EFTc"]
            f64 = evaluate(stmts, c, h0)
            mp = evaluate(stmts, c, h0, mp=mpmath)
            # replace rho_d (dust only) by the full (rho+P) of all non-DE species, mapped to an effective dust density
            rho_full = 2 * full / a ** 2

            def rho_over(ns, _r=rho_full):
                return mpmath.mpf(_r)

            mpf = evaluate(stmts, c, h0, mp=mpmath, override={"rho_d": rho_over})
            print(f"{a:5.2f} {c['EFTc'] / H ** 2:10.3e} {resid / dust:16.8f} {resid / full:10.6f} "
                  f"{c['EFT_mu1'] * norm:11.3e} {f64['C_EFT_mu1'] * norm:11.3e} {float(mp['C_EFT_mu1']) * norm:11.3e} "
                  f"{float(mp['C_EFT_mu2']) * norm:11.3e} {float(mpf['C_EFT_mu1']) * norm:12.3e} "
                  f"{float(mpf['C_EFT_mu2']) * norm:12.3e}")
        # c-scaling at a=0.5: scale c, cdot, cdotdot by s (keeping everything else) -> mu ~ s^-p ?
        c = [x for x in r["rows"] if abs(x["a"] - 0.5) < 1e-12][0]
        norm = (0.5 / c["adotoa"]) ** 2
        print("  c-scaling at a=0.5 (EFTc, EFTcdot, EFTcdotdot multiplied by s; mp60):")
        for s in (10.0, 1.0, 0.1, 0.01):
            cc = dict(c)
            for k in ("EFTc", "EFTcdot", "EFTcdotdot"):
                cc[k] = c[k] * s
            mp = evaluate(stmts, cc, h0, mp=mpmath)
            print(f"    s={s:6.2f}  mu1/H^2={float(mp['C_EFT_mu1']) * norm: .4e}  mu2/H^2={float(mp['C_EFT_mu2']) * norm: .4e}"
                  f"  3F2^2+F1F3 = {float(3 * mp['F2'] ** 2 + mp['F1'] * mp['F3']):.4e}  4c/a^2 = {4 * cc['EFTc'] / 0.25:.4e}")


if __name__ == "__main__":
    main()

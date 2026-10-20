"""Track t4a-param: numerical checks of every EFTCAMB 1D/2D parametrization.

Builds param_driver.f90 against the modules/library of a frozen audit build, instantiates every concrete
parametrization with several parameter sets through the abstract interface and compares

  * first..fourth derivatives with Richardson-extrapolated 5-point central differences of the next lower method
  * integral(a) with the EFTCAMB convention I(a) = a^2 exp(-3 int_1^a (1+w(a'))/a' da') (grhov_t = grhov*I(a),
    see 007p3_Designer_fR.f90:1100), computed by composite Gauss-Legendre in ln a from the driver's own value()
  * spline continuity of value/derivatives at the nodes and interpolation of the node values
  * parameter count / names / values round-trip, feedback, init_from_file keys, allocator flag -> type map

usage: python run_param_checks.py [--build A0|AD] [--out DIR] [--only substr]
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
HERE = Path(__file__).resolve().parent
CONDA_LIB = "/Users/marco/Software/miniforge3/lib"

DOC_ALLOC = {
    0: "zero_parametrization_1D",
    1: "constant_parametrization_1D",
    2: "linear_parametrization_1D",
    3: "power_law_parametrization_1D",
    4: "exponential_parametrization_1D",
    5: "taylorseries_parametrization_1D",
    6: "padeseries_parametrization_1D",
    7: "fourier_parametrization_1D",
    8: "steplog_parametrization_1D",
    9: "spline_parametrization_1D",
    10: "spline5_parametrization_1D",
    11: "exponential_parametrization_2_1D",
}


def build_driver(build: str, out: Path) -> Path:
    bdir = Path(build) if "/" in build else AUDIT / "build" / build
    exe = out / f"driver_{bdir.name}"
    src = HERE / "param_driver.f90"
    if exe.exists() and exe.stat().st_mtime > src.stat().st_mtime:
        return exe
    flags = ["-O0", "-g", "-cpp", "-ffree-line-length-none", "-fopenmp"]
    if bdir.name in ("AD", "ADD"):
        flags += ["-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
    cmd = [
        "gfortran",
        *flags,
        f"-I{bdir}/camb/Releaselib",
        f"-I{bdir}/forutils/Releaselib",
        str(src),
        "-o",
        str(exe),
        str(bdir / "camb/Releaselib/camblib.so"),
        f"-Wl,-rpath,{CONDA_LIB}",
    ]
    subprocess.run(cmd, check=True, cwd=out, capture_output=True)
    return exe


# ----------------------------------------------------------------------------------------------- cases


def A_POINTS(extra=()):
    base = [1e-8, 1e-6, 1e-4, 1e-2, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.999, 1.0]
    return sorted(set(base) | set(extra))


PHI_POINTS = [-2.0, -0.7, 0.3, 0.9, 1.0, 1.7, 2.9]


def cases() -> list[dict]:
    c: list[dict] = []

    def add(name, kind, key, fname, ini, pts, *, d3=True, d4=True, integ=False, nodes=(), npar=None, ptype=None,
            xvar="a", interior=None, flabel=None, int_conv="std", fb_noarg=True, neg_ok=False):
        c.append(
            dict(name=name, kind=kind, key=str(key), fname=fname, flabel=flabel or fname, ini=ini, pts=list(pts),
                 d3=d3, d4=d4, integ=integ, nodes=list(nodes), npar=npar, ptype=ptype, xvar=xvar,
                 interior=interior, int_conv=int_conv, fb_noarg=fb_noarg, neg_ok=neg_ok)
        )

    # neutral
    add("zero", "alloc", 0, "EFTOmega", {}, A_POINTS(), integ=True, npar=0, ptype=DOC_ALLOC[0])
    add("wlcdm", "type", "wlcdm", "EFTw", {}, A_POINTS(), integ=True, npar=0, ptype="wDE_LCDM_parametrization_1D")
    # constant
    for i, v in enumerate([0.3, -0.9, 0.0]):
        add(f"constant_{i}", "alloc", 1, "EFTw", {"EFTw0": v}, A_POINTS(), integ=True, npar=1, ptype=DOC_ALLOC[1])
    # linear
    for i, v in enumerate([0.5, -1.2]):
        add(f"linear_{i}", "alloc", 2, "EFTw", {"EFTw0": v}, A_POINTS(), integ=True, npar=1, ptype=DOC_ALLOC[2],
            neg_ok=True)
    # power law
    for i, (cc, n) in enumerate([(0.1, 2.0), (-0.3, 0.5), (0.2, 0.0), (0.05, -1.5), (-0.4, 3.7)]):
        add(f"powerlaw_{i}", "alloc", 3, "EFTw", {"EFTw0": cc, "EFTwExp": n}, A_POINTS(), integ=True, npar=2,
            ptype=DOC_ALLOC[3])
    # exponential exp(c a^n)-1 (integral not implemented -> MpiStop, tested separately)
    for i, (cc, n) in enumerate([(0.1, 1.0), (-0.5, 2.0), (0.3, 0.5)]):
        add(f"exponential_{i}", "alloc", 4, "EFTOmega", {"EFTOmega0": cc, "EFTOmegaExp": n}, A_POINTS(), npar=2,
            ptype=DOC_ALLOC[4])
    # Taylor series (allocator 5)
    for i, (order, a0, coeffs) in enumerate([(0, 0.0, [0.3]), (3, 0.5, [0.1, -0.2, 0.3, 0.4]),
                                             (6, 1.0, [0.01, 0.02, -0.03, 0.04, 0.05, -0.06, 0.07])]):
        ini = {"EFTOmega_Taylor_order": order, "EFTOmegaa0": a0}
        ini.update({f"EFTOmega{j}": v for j, v in enumerate(coeffs)})
        add(f"taylorseries_{i}", "alloc", 5, "EFTOmega", ini, A_POINTS(), npar=order + 2, ptype=DOC_ALLOC[5],
            neg_ok=True)
    # Pade
    for i, (nn, nd, a0, up, down) in enumerate([
        (2, 1, 1.0, [0.1, 0.2, -0.3], [0.4]),
        (1, 2, 0.3, [0.05, 0.1], [0.2, -0.3]),
        (3, 3, 0.0, [0.1, -0.1, 0.2, 0.05], [0.3, 0.1, -0.2]),
        (4, 4, 0.5, [0.1, -0.1, 0.2, 0.05, 0.3], [0.3, 0.1, -0.2, 0.25]),
    ]):
        ini = {"EFTOmega_Pade_order_N": nn, "EFTOmega_Pade_order_D": nd, "EFTOmegaa0": a0}
        ini.update({f"EFTOmegaN{j}": v for j, v in enumerate(up)})
        ini.update({f"EFTOmegaD{j + 1}": v for j, v in enumerate(down)})
        add(f"pade_{i}", "alloc", 6, "EFTOmega", ini, A_POINTS(), npar=nn + nd + 2, ptype=DOC_ALLOC[6], neg_ok=True)
    # Fourier
    for i, (nc, ns, mina0, a0, ph, ca, sb) in enumerate([
        (2, 1, False, 0.1, 0.1, [0.3, -0.2], [0.15]),
        (2, 2, True, 0.05, 0.0, [0.3, -0.2], [0.15, 0.1]),
        (1, 0, False, 0.2, 0.25, [0.5], []),
        (0, 0, True, 0.2, 0.0, [], []),
    ]):
        ini = {"EFTOmega_Fourier_order_cos": nc, "EFTOmega_Fourier_order_sin": ns,
               "EFTOmega_min_a0": "T" if mina0 else "F", "EFTOmega_a0": a0, "EFTOmega_phase": ph}
        ini.update({f"EFTOmega_a{j + 1}": v for j, v in enumerate(ca)})
        ini.update({f"EFTOmega_b{j + 1}": v for j, v in enumerate(sb)})
        add(f"fourier_{i}", "alloc", 7, "EFTOmega", ini, A_POINTS(), npar=2 + nc + ns, ptype=DOC_ALLOC[7],
            neg_ok=True)
    # steplog (wDE-like, integral implemented)
    for i, (v1, v2, at, dl) in enumerate([(-1.0, -0.8, 0.3, 0.5), (-0.9, -1.1, 0.5, 0.1), (-0.95, -0.95, 0.2, 1.0)]):
        ini = {"EFTw_v1": v1, "EFTw_v2": v2, "EFTw_at": at, "EFTw_delta": dl}
        extra = [at, at * math.exp(dl), at * math.exp(-dl)]
        add(f"steplog_{i}", "alloc", 8, "EFTw", ini, A_POINTS(extra), integ=True, npar=4, ptype=DOC_ALLOC[8])
    # cubic spline (allocator 9)
    xs6 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    xs7 = [0.01, 0.1, 0.25, 0.4, 0.6, 0.8, 1.0]
    for i, (xs, vs, null) in enumerate([
        (xs6, [0.1 + 0.5 * x for x in xs6], 0.0),  # linear data -> exact spline is linear
        (xs7, [-1.0 + 0.3 * x * x for x in xs7], -1.0),  # w-like data, null value -1
        (xs7, [0.2 * math.sin(4 * x) for x in xs7], 0.0),
    ]):
        n = len(xs)
        ini = {"EFTw_Spline_Pixels": n, "EFTw_null_value": null}
        ini.update({f"EFTwx{j + 1}": x for j, x in enumerate(xs)})
        ini.update({f"EFTwv{j + 1}": v for j, v in enumerate(vs)})
        interior = [0.5 * (xs[j] + xs[j + 1]) for j in range(n - 1)] + [xs[0] + 0.3 * (xs[1] - xs[0]),
                                                                       xs[-1] - 0.3 * (xs[-1] - xs[-2])]
        add(f"spline_{i}", "alloc", 9, "EFTw", ini, sorted(interior), integ=True, nodes=xs, npar=2 * n,
            ptype=DOC_ALLOC[9], interior=(xs, vs, null))
    # spline5 (allocator 10)
    xs8 = [1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, 0.6, 0.8, 1.0]
    for i, (order, logx, logy, vs, null) in enumerate([
        (5, False, False, [0.2 * math.sin(4 * x) + 0.1 for x in xs8], 0.0),
        (5, True, False, [-1.0 + 0.1 * x for x in xs8], -1.0),  # w-like, needed for integral
        (3, False, False, [0.2 * math.sin(4 * x) + 0.1 for x in xs8], 0.0),
        (5, True, True, [0.5 + 0.3 * x for x in xs8], 0.0),
        (5, False, False, [0.1 + 0.5 * x for x in xs8], 0.0),  # linear data -> exact
    ]):
        n = len(xs8)
        ini = {"EFTw_Spline_Order": order, "EFTw_Logx_Flag": "T" if logx else "F",
               "EFTw_Logy_Flag": "T" if logy else "F", "EFTw_Spline_Pixels": n, "EFTw_null_value": null}
        ini.update({f"EFTwx{j + 1}": x for j, x in enumerate(xs8)})
        ini.update({f"EFTwv{j + 1}": v for j, v in enumerate(vs)})
        if logx:
            interior = [math.sqrt(xs8[j] * xs8[j + 1]) for j in range(n - 1)]
        else:
            interior = [0.5 * (xs8[j] + xs8[j + 1]) for j in range(n - 1)]
        add(f"spline5_{i}", "alloc", 10, "EFTw", ini, sorted(interior), d3=(order == 5), d4=(order == 5),
            integ=(logx and not logy), nodes=xs8, npar=2 * n, ptype=DOC_ALLOC[10], interior=(xs8, vs, null))
    # exp2 c*exp(-l x): quintessence potential variable phi
    for i, (cc, ll) in enumerate([(0.7, 2.0), (1e-3, -0.5)]):
        add(f"exp2_{i}", "alloc", 11, "V", {"V0": cc, "V1": ll}, PHI_POINTS + [1e-8], npar=2,
            ptype=DOC_ALLOC[11], xvar="phi")
    # wDE-only types allocated by name in the designer models
    for i, (w0, wa) in enumerate([(-0.9, 0.3), (-1.1, -0.5)]):
        add(f"cpl_{i}", "type", "cpl", "EFTw", {"EFTw0": w0, "EFTw1": wa}, A_POINTS(), integ=True, npar=2,
            ptype="CPL_parametrization_1D", neg_ok=True)
    for i, (w0, wa, wn) in enumerate([(-0.9, 0.3, 2.0), (-1.0, 0.2, 1.0), (-0.95, 0.4, 3.5), (-0.9, 0.3, 1.0 + 1e-9),
                                      (-0.9, 0.3, 0.5), (-0.9, 0.3, 1.0 + 1e-6)]):
        add(f"jbp_{i}", "type", "jbp", "EFTw", {"EFTw0": w0, "EFTw1": wa, "EFTw2": wn}, A_POINTS(), integ=True,
            npar=3, ptype="JBP_parametrization_1D")
    for i, (w0, wa, wat) in enumerate([(-0.9, 0.5, 0.4), (-1.0, -0.3, 1.2)]):
        add(f"turning_{i}", "type", "turning", "EFTw", {"EFTw0": w0, "EFTw1": wa, "EFTw2": wat},
            A_POINTS([wat] if wat < 1 else []), integ=True, npar=3, ptype="turning_point_parametrization_1D",
            neg_ok=True)
    add("taylor_0", "type", "taylor", "EFTw", {"EFTw0": -0.9, "EFTw1": 0.2, "EFTw2": 0.3, "EFTw3": -0.4},
        A_POINTS(), integ=True, npar=4, ptype="taylor_parametrization_1D", neg_ok=True)
    # quintessence potentials (variable phi, no integral)
    add("dexp_0", "type", "dexp", "V", {"V0": 0.7, "V1": 2.0, "V2": 0.3, "V3": -1.0}, PHI_POINTS + [1e-8], npar=4,
        ptype="double_exponential_parametrization_1D", xvar="phi")
    add("cosine_0", "type", "cosine", "V", {"V0": 0.7, "V1": 3.0, "V2": 0.4, "V3": 1.2}, PHI_POINTS + [1e-8],
        npar=4, ptype="cosine_parametrization_1D", xvar="phi")
    for i, (f0, al, n) in enumerate([(0.7, 2.0, 1.0), (0.7, 2.0, 2.5), (0.7, 1.3, 0.5), (0.7, 2.0, 3.0)]):
        add(f"axion_{i}", "type", "axion", "V", {"V0": f0, "V1": al, "V2": n}, [0.3, 0.9, 1.0, 1.7, 2.9, -0.7],
            npar=3, ptype="axion_parametrization_1D", xvar="phi")
    add("powsum_0", "type", "powsum", "V",
        {"V0": 0.1, "V1": 1.0, "V2": 2.0, "V3": -0.5, "V4": 3.0, "V5": 0.2, "V6": -1.0},
        [0.05, 0.3, 0.9, 1.0, 1.7, 2.9], npar=7, ptype="power_law_sum_parametrization_1D", xvar="phi")
    # step (not reachable from the allocator; integral implemented)
    for i, (v1, v2, at, dl) in enumerate([(-1.0, -0.8, 0.3, 0.1), (-0.9, -1.1, 0.5, 0.5), (-1.0, -0.9, 0.4, 1.0)]):
        add(f"step_{i}", "type", "step", "EFTw", {"EFTw_v1": v1, "EFTw_v2": v2, "EFTw_at": at, "EFTw_delta": dl},
            A_POINTS([at, at + dl, max(at - dl, 0.02)]), integ=True, npar=4, ptype="step_parametrization_1D",
            neg_ok=True)
    return c


# ----------------------------------------------------------------------------------------------- running


def write_ini(path: Path, d: dict):
    lines = [f"{k} = {v}" for k, v in d.items()]
    path.write_text("\n".join(lines) + "\n")


H_REL = 2e-3


def stencil(x0: float, case: dict) -> tuple[float, list[float]]:
    if case["xvar"] == "phi" or case["neg_ok"]:
        h = H_REL * max(abs(x0), 0.05)
    else:
        h = H_REL * x0
    if case["nodes"]:
        nodes = np.array(case["nodes"], dtype=float)
        if case["key"] == "10" and case["ini"].get("EFTw_Logx_Flag") == "T":
            pass
        d = np.min(np.abs(nodes - x0))
        h = min(h, d / 6.0)
    pts = [x0 + m * h for m in range(-4, 5)] + [x0 + m * h / 2 for m in (-3, -1, 1, 3)]
    return h, pts


def gl_grid(lo: float, hi: float, breaks: list[float], panel: float = 0.02, order: int = 12):
    """Composite Gauss-Legendre nodes on [lo, hi] (in N = ln a) with panel edges at breaks."""
    edges = sorted(set([lo, hi] + [b for b in breaks if lo < b < hi]))
    xg, wg = np.polynomial.legendre.leggauss(order)
    nodes, weights, panel_of = [], [], []
    pe = []
    for e0, e1 in zip(edges[:-1], edges[1:]):
        m = max(1, int(math.ceil((e1 - e0) / panel)))
        sub = np.linspace(e0, e1, m + 1)
        for s0, s1 in zip(sub[:-1], sub[1:]):
            pe.append((s0, s1))
            nodes.extend(0.5 * (s1 - s0) * xg + 0.5 * (s1 + s0))
            weights.extend(0.5 * (s1 - s0) * wg)
    return np.array(nodes), np.array(weights), pe


def run_case(exe: Path, case: dict, wdir: Path) -> dict:
    cdir = wdir / "cases"
    cdir.mkdir(parents=True, exist_ok=True)
    ini_path = cdir / f"{case['name']}.ini"
    write_ini(ini_path, case["ini"])
    allpts: list[float] = []
    tags: list[tuple] = []
    for x0 in case["pts"]:
        h, pts = stencil(x0, case)
        tags.append(("stencil", x0, h, len(allpts)))
        allpts.extend(pts)
    # node continuity points
    for xn in case["nodes"]:
        for eps in (-1e-9, 0.0, 1e-9):
            allpts.append(xn + eps * max(abs(xn), 1e-3))
    # integral grid in N = ln a
    int_info = None
    if case["integ"]:
        test_a = [a for a in case["pts"] if 0 < a <= 1.0]
        breaks = [math.log(a) for a in test_a]
        if case["nodes"]:
            breaks += [math.log(x) for x in case["nodes"] if x > 0]
        nmin = math.log(min(test_a))
        nodes_n, w_n, _ = gl_grid(nmin, 0.0, breaks)
        int_info = (len(allpts), nodes_n, w_n, test_a)
        allpts.extend(np.exp(nodes_n).tolist())
    case_file = cdir / f"{case['name']}.case"
    with open(case_file, "w") as f:
        f.write(f"{case['kind']}\n{case['key']}\n{case['fname']}\n{case['flabel']}\n{ini_path}\n")
        tf = lambda b: "T" if b else "F"  # noqa: E731
        f.write(f"{tf(case['d3'])} {tf(case['d4'])} {tf(case['integ'])} {tf(case['fb_noarg'])}\n")
        f.write(f"{len(allpts)}\n")
        for x in allpts:
            f.write(f"{x:.17e}\n")
    p = subprocess.run([str(exe), str(case_file)], capture_output=True, text=True, timeout=600)
    (cdir / f"{case['name']}.out").write_text(p.stdout + "\nSTDERR\n" + p.stderr)
    res = parse(p.stdout)
    res["rc"] = p.returncode
    res["stderr_tail"] = p.stderr[-2000:]
    res["tags"] = tags
    res["int_info"] = int_info
    res["nodes_start"] = len(case["pts"]) and tags[-1][3] + len(stencil(case["pts"][-1], case)[1])
    return res


def parse(out: str) -> dict:
    r: dict = {"params": [], "data": [], "feedback": [], "type": None, "tail": out.splitlines()[-3:]}
    infb = False
    for line in out.splitlines():
        if line.startswith("FEEDBACK_BEGIN"):
            infb = True
            continue
        if line.startswith("FEEDBACK_END"):
            infb = False
            continue
        if infb:
            r["feedback"].append(line)
            continue
        tok = line.split()
        if not tok:
            continue
        if tok[0] == "TYPE":
            r["type"] = tok[1]
        elif tok[0] in ("ALLOC_ERROR", "INITFUNC_ERROR", "INIT_ERROR"):
            r[tok[0]] = int(tok[1])
        elif tok[0] == "NPARAM":
            r["nparam"] = (int(tok[1]), int(tok[2]))
        elif tok[0] == "PARAM":
            rest = line.split("|")
            r["params"].append(dict(i=int(tok[1]), value=float(tok[2]), name=rest[1], short=rest[2], latex=rest[3]))
        elif tok[0] == "ROUNDTRIP":
            r["roundtrip"] = (float(tok[1]), float(tok[2]))
        elif tok[0] == "D":
            r["data"].append([float(t) for t in tok[1:8]] + [tok[8] == "T"])
        elif tok[0] == "DONE_EVAL":
            r["done"] = True
        elif tok[0] == "DONE":
            r["done_all"] = True
    return r


# ----------------------------------------------------------------------------------------------- checks


def fd_deriv(fvals: np.ndarray, h: float) -> tuple[float, float]:
    """fvals at x0 + m h, m=-4..4 then x0 + m h/2, m=-3,-1,1,3. Returns Richardson value and error estimate."""
    f = {m: fvals[m + 4] for m in range(-4, 5)}
    g = {-3: fvals[9], -1: fvals[10], 1: fvals[11], 3: fvals[12]}
    # 5-point with step 2h, h and h/2  (O(s^4))
    d_2h = (8 * (f[2] - f[-2]) - (f[4] - f[-4])) / (24 * h)
    d_h = (8 * (f[1] - f[-1]) - (f[2] - f[-2])) / (12 * h)
    d_h2 = (8 * (g[1] - g[-1]) - (f[1] - f[-1])) / (6 * h)
    r1 = (16 * d_h - d_2h) / 15
    r2 = (16 * d_h2 - d_h) / 15
    r = (64 * r2 - r1) / 63
    return r, abs(r - r2) + abs(r2 - r1) * 0.1


METHODS = ["value", "d1", "d2", "d3", "d4"]


def analyse(case: dict, res: dict) -> dict:
    out: dict = {"name": case["name"], "status": "ok", "issues": []}
    if not res.get("done"):
        out["status"] = f"rc={res['rc']}"
        out["issues"].append(f"driver did not finish: rc={res['rc']} {res['stderr_tail'][-400:]!r}")
        return out
    data = np.array([row[:7] for row in res["data"]], dtype=float)
    # metadata checks
    if case["ptype"] and res["type"] != case["ptype"]:
        out["issues"].append(f"type {res['type']} != expected {case['ptype']}")
    if case["npar"] is not None and res.get("nparam", (None,))[0] != case["npar"]:
        out["issues"].append(f"parameter_number {res.get('nparam')} != expected {case['npar']}")
    if len(res["params"]) != res.get("nparam", (0,))[0]:
        out["issues"].append("number of PARAM lines != parameter_number")
    names = [p["name"] for p in res["params"]]
    if len(set(names)) != len(names):
        out["issues"].append(f"duplicate parameter names {names}")
    for p in res["params"]:
        if p["name"] != p["short"]:
            out["issues"].append(f"parameter name truncated to EFT_names_max_length: {p['name']} -> {p['short']}")
        key = p["name"]
        if key in case["ini"]:
            want = case["ini"][key]
            if isinstance(want, (int, float)) and not isinstance(want, bool):
                # spline5 logx/logy store log(x) / log(v): compare via transformation
                got = p["value"]
                if case["key"] == "10" and case["ini"].get("EFTw_Logx_Flag") == "T" and key.startswith("EFTwx"):
                    want = math.log(want)
                if case["key"] == "10" and case["ini"].get("EFTw_Logy_Flag") == "T" and key.startswith("EFTwv"):
                    want = math.log(want)
                if case["key"] == "7" and case["ini"].get("EFTOmega_min_a0") == "T" and key == "EFTOmega_a0":
                    pass
                if abs(got - want) > 1e-12 * max(1, abs(want)):
                    out["issues"].append(f"parameter_value({key}) = {got} but ini has {want}")
        else:
            out["issues"].append(f"parameter key {key!r} not set in test ini (name mismatch with docs?)")
    rt = res.get("roundtrip")
    if not res.get("done_all"):
        out["issues"].append("init_parameters(parameter_value(:)) round trip aborted the program: "
                             + res["stderr_tail"][-200:].replace("\n", " ") + " | stdout tail: "
                             + " ".join(res.get("tail", [])))
    if rt and not (rt[0] == rt[1] or (math.isnan(rt[0]) and math.isnan(rt[1]))):
        out["issues"].append(f"init_parameters(parameter_value(:)) round trip changes value: {rt}")
    fb = "\n".join(res["feedback"])
    out["feedback_lines"] = len(res["feedback"])

    # derivative checks
    worst = {m: (0.0, None) for m in METHODS[1:]}
    nonfinite = []
    have = {"d1": True, "d2": True, "d3": case["d3"], "d4": case["d4"]}
    per_point = []
    for tag in res["tags"]:
        _, x0, h, i0 = tag
        blk = data[i0:i0 + 13]
        row = data[i0 + 4]
        if not np.all(np.isfinite(row[1:6])):
            nonfinite.append((x0, row[1:6].tolist()))
        pp = {"x0": x0}
        for k in range(1, 5):
            m = METHODS[k]
            if not have[m]:
                continue
            ref, est = fd_deriv(blk[:, k], h)  # derivative of method k-1 (column k is value..d3)
            code = row[k + 1]
            pp[m] = (code, ref, est)
            # roundoff budget of the Richardson stencil (~3.4/h times the evaluation error, taken as 8 ulp)
            fmax = np.max(np.abs(blk[:, k])) if np.all(np.isfinite(blk[:, k])) else 0.0
            pp["noise_" + m] = 3.4 * 8 * 2.2e-16 * fmax / h
        per_point.append(pp)
    # normalise by typical magnitude of each method over the points with x0 in [0.01, 1] (a) / all (phi)
    for m in METHODS[1:]:
        if not have[m]:
            continue
        refs = [abs(pp[m][1]) for pp in per_point if m in pp and np.isfinite(pp[m][1])
                and (case["xvar"] == "phi" or 0.01 <= pp["x0"] <= 1.0)]
        floor = 1e-3 * max(refs) if refs else 0.0
        vmax = max([abs(r[1]) for r in data[:, :2] if np.isfinite(r[1])] or [0.0])
        floor = max(floor, 1e-8 * vmax, 1e-12)
        for pp in per_point:
            if m not in pp:
                continue
            code, ref, est = pp[m]
            if not (np.isfinite(code) and np.isfinite(ref)):
                if np.isfinite(ref) != np.isfinite(code):
                    worst[m] = (math.inf, pp["x0"])
                continue
            err = max(abs(code - ref) - 3 * est - pp["noise_" + m], 0.0) / max(abs(ref), floor)
            if err > worst[m][0]:
                worst[m] = (err, pp["x0"], code, ref)
    out["worst"] = {m: worst[m] for m in METHODS[1:] if have[m]}
    out["nonfinite"] = nonfinite
    tol = {"d1": 1e-6, "d2": 1e-4, "d3": 1e-4, "d4": 1e-4}
    for m, w in out["worst"].items():
        if w[0] > tol[m]:
            out["issues"].append(f"{m}: worst rel err {w[0]:.3e} at x={w[1]} (code {w[2] if len(w) > 2 else None}, "
                                 f"fd {w[3] if len(w) > 3 else None})")

    # node checks (splines)
    if case["nodes"]:
        start = res["tags"][-1][3] + 13
        xs, vs, null = case["interior"]
        cont = {m: 0.0 for m in METHODS}
        interp = 0.0
        node_report = []
        for j, xn in enumerate(case["nodes"]):
            lo, mid, hi = data[start + 3 * j], data[start + 3 * j + 1], data[start + 3 * j + 2]
            interp_err = abs(mid[1] - vs[j]) / max(abs(vs[j]), 1e-3)
            node_report.append((xn, mid[1], vs[j]))
            interp = max(interp, interp_err)
            if 0 < j < len(case["nodes"]) - 1:
                for k, m in enumerate(METHODS):
                    if m in ("d3", "d4") and (not have[m] or case["key"] == "9" or
                                              (case["key"] == "10" and case["ini"].get("EFTw_Spline_Order") == 3)):
                        continue
                    col = data[:, k + 1]
                    cmax = np.max(np.abs(col[np.isfinite(col)])) if np.any(np.isfinite(col)) else 0.0
                    vmax = np.max(np.abs(data[:, 1][np.isfinite(data[:, 1])]))
                    scale = max(abs(lo[k + 1]), abs(hi[k + 1]), 1e-4 * cmax, 1e-8 * vmax, 1e-10)
                    cont[m] = max(cont[m], abs(hi[k + 1] - lo[k + 1]) / scale)
        out["node_interp_err"] = interp
        out["node_values"] = node_report
        out["node_cont"] = cont
        if interp > 1e-10:
            out["issues"].append(f"spline does not reproduce node values: worst rel {interp:.3e}; nodes "
                                 f"(x, value(x), v): {[(a, round(b, 12), c) for a, b, c in node_report]}")
        for m, v in cont.items():
            lim = 1e-6 if m in ("value", "d1", "d2") else 1e-6
            if v > lim:
                out["issues"].append(f"discontinuity of {m} at interior nodes: max rel jump {v:.3e}")
    # integral check
    if case["integ"]:
        i0, nodes_n, w_n, test_a = res["int_info"]
        wvals = data[i0:i0 + len(nodes_n), 1]
        ints = []
        worst_i = (0.0, None)
        for a in test_a:
            # value at the plain point x0 = a is the centre of its stencil
            tag = next(t for t in res["tags"] if t[1] == a)
            code = data[tag[3] + 4, 6]
            na = math.log(a)
            sel = nodes_n >= na - 1e-15
            integral = np.sum(w_n[sel] * (1.0 + wvals[sel]))  # int_{ln a}^0 (1+w) dN
            logref = 2.0 * math.log(a) + 3.0 * integral
            ref = math.exp(logref) if logref < 709 else math.inf
            if code > 0 and np.isfinite(code):
                err = abs(math.expm1(math.log(code) - logref))
            elif (code == math.inf and logref > 709) or (code == 0 and logref < -745):
                err = 0.0
            else:
                err = math.inf
            ints.append((a, code, ref, err))
            if err > worst_i[0] or (math.isnan(err)):
                worst_i = (err, a, code, ref)
        out["integral"] = ints
        out["worst_integral"] = worst_i
        if not (worst_i[0] <= 1e-6):
            out["issues"].append(f"integral: worst rel err {worst_i[0]:.3e} at a={worst_i[1]} (code {worst_i[2]}, "
                                 f"ref {worst_i[3]})")
    if nonfinite:
        out["issues"].append(f"non-finite outputs at {[(x, v) for x, v in nonfinite][:4]}")
    if out["issues"]:
        out["status"] = "issues"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default="A0")
    ap.add_argument("--out", default=str(AUDIT / "work" / "t4a-param"))
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    wdir = Path(args.out) / Path(args.build).name
    wdir.mkdir(parents=True, exist_ok=True)
    exe = build_driver(args.build, Path(args.out))
    results = []
    for case in cases():
        if args.only and args.only not in case["name"]:
            continue
        res = run_case(exe, case, wdir)
        a = analyse(case, res)
        results.append(a)
        flag = "OK " if a["status"] == "ok" else "BAD"
        worst = a.get("worst", {})
        ws = " ".join(f"{m}={v[0]:.1e}" for m, v in worst.items())
        wi = a.get("worst_integral")
        if wi:
            ws += f" int={wi[0]:.1e}"
        print(f"{flag} {case['name']:16s} {ws}")
        for iss in a["issues"]:
            print(f"      - {iss[:400]}")
    (wdir / "results.json").write_text(json.dumps(results, indent=1, default=str))


if __name__ == "__main__":
    sys.exit(main())

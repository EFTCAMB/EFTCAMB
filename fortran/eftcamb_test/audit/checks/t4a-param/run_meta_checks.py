"""Track t4a-param: allocator map, parameter-name length, 2D parametrizations, feedback, Fourier minimum.

usage: python run_meta_checks.py [--build A0|AD]
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_param_checks as rp  # noqa: E402


def run(exe: Path, wdir: Path, name: str, kind: str, key: str, fname: str, ini: dict, pts, flags=(True, True,
        False, True), pts2=None):
    cdir = wdir / "meta"
    cdir.mkdir(parents=True, exist_ok=True)
    ini_path = cdir / f"{name}.ini"
    rp.write_ini(ini_path, ini)
    cf = cdir / f"{name}.case"
    tf = lambda b: "T" if b else "F"  # noqa: E731
    with open(cf, "w") as f:
        f.write(f"{kind}\n{key}\n{fname}\n{fname}\n{ini_path}\n{' '.join(tf(b) for b in flags)}\n{len(pts)}\n")
        for i, x in enumerate(pts):
            f.write(f"{x:.17e} {pts2[i]:.17e}\n" if pts2 is not None else f"{x:.17e}\n")
    p = subprocess.run([str(exe), str(cf)], capture_output=True, text=True, timeout=300)
    (cdir / f"{name}.out").write_text(p.stdout + "\nSTDERR\n" + p.stderr)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default="A0")
    args = ap.parse_args()
    out = rp.AUDIT / "work" / "t4a-param"
    exe = rp.build_driver(args.build, out)
    wdir = out / Path(args.build).name
    report: dict = {}

    # 1) allocator: flag -> dynamic type, invalid flags -> eft_error = 1
    amap = {}
    for flag in list(range(0, 12)) + [12, -1, 99]:
        p = run(exe, wdir, f"alloc_{flag}", "alloc", str(flag), "EFTOmega", {}, [0.5], flags=(False, False, False,
                                                                                               True))
        r = rp.parse(p.stdout)
        amap[flag] = dict(type=r["type"], alloc_error=r.get("ALLOC_ERROR"), rc=p.returncode,
                          nparam=r.get("nparam"), ok=(r["type"] == rp.DOC_ALLOC.get(flag)) if flag in rp.DOC_ALLOC
                          else (r.get("ALLOC_ERROR") == 1))
    report["allocator"] = amap
    print("allocator map:")
    for k, v in amap.items():
        print(f"   {k:3d} -> {v}")

    # 2) EFT_names_max_length = 20: long function names used by the 3+1 designer
    trunc = {}
    for fname in ["Designer3p1sigma1_ODE", "Designer3p1rho1N_ODE", "Designer3p1rho2_ODE", "RPHkineticity_ODE"]:
        for flag, keys in [(3, {f"{fname}0": 0.7, f"{fname}Exp": 2.0}),
                           (8, {f"{fname}_v1": 0.1, f"{fname}_v2": 0.2, f"{fname}_at": 0.3, f"{fname}_delta": 0.4})]:
            p = run(exe, wdir, f"trunc_{fname}_{flag}", "alloc", str(flag), fname, keys, [0.5],
                    flags=(False, False, False, True))
            r = rp.parse(p.stdout)
            got = {pp["name"]: (pp["value"], pp["short"]) for pp in r["params"]}
            want = {k: v for k, v in keys.items()}
            bad = {k: dict(read=got.get(k, (None, None))[0], key_used=got.get(k, (None, None))[1], ini=v)
                   for k, v in want.items() if got.get(k, (None,))[0] != v}
            trunc[f"{fname}/flag{flag}"] = bad
    report["name_truncation"] = trunc
    print("name truncation (keys not read back):")
    for k, v in trunc.items():
        print(f"   {k}: {v}")

    # 3) 2D parametrizations
    xs = np.array([1e-8, 0.1, 0.5, 1.0, 2.0, -0.3])
    ys = np.array([0.3, -1.0, 0.7, 1.0, 1e-8, 2.0])
    two = {}
    for key, ini, f in [("zero", {}, lambda x, y: (0, 0, 0, 0, 0, 0)),
                        ("constant", {"F0": 0.37}, lambda x, y: (0.37, 0, 0, 0, 0, 0)),
                        ("bilinear", {"F0": 0.5, "F1": -1.5}, lambda x, y: (0.5 * x - 1.5 * y, 0.5, -1.5, 0, 0, 0))]:
        p = run(exe, wdir, f"2d_{key}", "2d", key, "F", ini, xs, flags=(False, False, False, False), pts2=ys)
        rows = [line.split()[1:] for line in p.stdout.splitlines() if line.startswith("D2 ")]
        err = 0.0
        for row in rows:
            x, y, *vals = map(float, row)
            ref = f(x, y)
            err = max(err, max(abs(a - b) for a, b in zip(vals, ref)))
        r = rp.parse(p.stdout)
        params = [(pp["name"], pp["value"]) for pp in r["params"]]
        p2 = run(exe, wdir, f"2d_{key}_fbnoarg", "2d", key, "F", ini, xs, flags=(False, False, False, True), pts2=ys)
        two[key] = dict(max_abs_err=err, rc=p.returncode, nparam=r.get("nparam"), params=params,
                        roundtrip=r.get("roundtrip"), feedback_noarg_rc=p2.returncode,
                        feedback_noarg_stderr=p2.stderr[-300:])
    report["2d"] = two
    print("2D:")
    for k, v in two.items():
        print(f"   {k}: {v}")

    # 4) Fourier min_a0: minimum over a period must equal a0
    fm = {}
    for i, (nc, ns, a0, ph, ca, sb) in enumerate([(2, 2, 0.05, 0.0, [0.3, -0.2], [0.15, 0.1]),
                                                  (3, 1, -0.2, 0.37, [0.1, 0.4, -0.2], [0.3]),
                                                  (0, 0, 0.2, 0.0, [], [])]):
        ini = {"EFTOmega_Fourier_order_cos": nc, "EFTOmega_Fourier_order_sin": ns, "EFTOmega_min_a0": "T",
               "EFTOmega_a0": a0, "EFTOmega_phase": ph}
        ini.update({f"EFTOmega_a{j + 1}": v for j, v in enumerate(ca)})
        ini.update({f"EFTOmega_b{j + 1}": v for j, v in enumerate(sb)})
        grid = ph + np.linspace(0, 1, 20001)
        p = run(exe, wdir, f"fourier_min_{i}", "alloc", "7", "EFTOmega", ini, grid, flags=(False, False, False, True))
        vals = np.array([float(line.split()[2]) for line in p.stdout.splitlines() if line.startswith("D ")])
        fm[i] = dict(a0=a0, grid_min=float(vals.min()), err=float(vals.min() - a0))
    report["fourier_min"] = fm
    print("Fourier min_a0:", fm)

    # 5) integral of functions without an analytic integral, and MpiStop paths
    noint = {}
    for name, kind, key, ini in [
        ("exponential", "alloc", "4", {"EFTw0": 0.1, "EFTwExp": 1.0}),
        ("taylorseries", "alloc", "5", {"EFTw_Taylor_order": 1, "EFTw0": -0.9}),
        ("pade", "alloc", "6", {"EFTwN0": -0.9}),
        ("fourier", "alloc", "7", {"EFTw_a0": -0.9}),
        ("exp2", "alloc", "11", {"EFTw0": -0.9}),
        ("dexp", "type", "dexp", {"EFTw0": -0.9}),
        ("cosine", "type", "cosine", {"EFTw0": -0.9}),
        ("axion", "type", "axion", {"EFTw0": -0.9}),
        ("powsum", "type", "powsum", {"EFTw0": -0.9}),
    ]:
        p = run(exe, wdir, f"noint_{name}", kind, key, "EFTw", ini, [0.5, 1.0], flags=(True, True, True, True))
        rows = [line.split() for line in p.stdout.splitlines() if line.startswith("D ")]
        noint[name] = dict(rc=p.returncode, integral=[float(r[7]) for r in rows],
                           stopped=("not implemented" in p.stdout), tail=p.stdout.splitlines()[-2:])
    report["integral_unimplemented"] = noint
    print("integral() without analytic solution:")
    for k, v in noint.items():
        print(f"   {k}: {v}")
    (wdir / "meta_results.json").write_text(json.dumps(report, indent=1, default=str))


if __name__ == "__main__":
    main()

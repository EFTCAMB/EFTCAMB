"""Verifier: search for a stability-verdict flip caused by the tripled gpinu_tot (buggy tree vs fixed tree).

For (model ini, parameter key, value grid, stability-flag config): evaluates the end-to-end stability verdict
(camb.get_background with EFTCAMB_skip_stability = F) on the grid in both trees, then bisects every verdict change in
each tree separately to locate the boundaries, and reports where the two trees disagree.

Usage: python gpinu_param_scan.py OUT.json TREE_BUGGY TREE_FIXED
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

PAR = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")
AUD = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/inis")

FLAGS = {
    "ini": {},
    "mass_math": {"EFT_mass_math_stability": True},
    "mass": {"EFT_mass_stability": True},
    "positivity": {"EFT_positivity_bounds": True},
    "minkowski": {"EFT_minkowski_limit": True},
}

WORKER = r"""
import json, sys, tempfile, os
sys.path.insert(0, sys.argv[1])
import camb
base, key, vals, flags = sys.argv[2], sys.argv[3], json.loads(sys.argv[4]), json.loads(sys.argv[5])
out = []
d = tempfile.mkdtemp()
for v in vals:
    ini = os.path.join(d, "m.ini")
    with open(ini, "w") as f:
        f.write(f"DEFAULT({base})\n{key} = {v!r}\noutput_root = {d}/out\nfeedback_level = 0\nEFTCAMB_write_background = F\n")
    cp = camb.read_ini(ini, no_validate=True)
    cp.EFTCAMB.EFTCAMB_skip_stability = False
    cp.EFTCAMB.EFTCAMB_feedback_level = 0
    for k, x in flags.items():
        setattr(cp.EFTCAMB, k, x)
    try:
        camb.get_background(cp)
        out.append("S")
    except Exception as e:  # noqa
        msg = str(e).strip().splitlines()[-1]
        out.append("U" if "unstable" in msg else "E")
print("RESULT " + json.dumps(out))
"""


def verdicts(tree, base, key, vals, flags):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with tempfile.TemporaryDirectory() as d, harness.slot(1):
        p = subprocess.run([sys.executable, "-c", WORKER, tree, str(base), key, json.dumps([float(v) for v in vals]),
                            json.dumps(flags)], capture_output=True, text=True, cwd=d, env=env,
                           preexec_fn=lambda: os.nice(10), timeout=900)
    res = [json.loads(line[7:]) for line in p.stdout.splitlines() if line.startswith("RESULT ")]
    if not res:
        raise RuntimeError(p.stdout[-800:] + p.stderr[-800:])
    return res[0]


def bisect(tree, base, key, lo, hi, vlo, flags, nsteps=int(os.environ.get("NBISECT", "24"))):
    for _ in range(nsteps):
        mid = 0.5 * (lo + hi)
        vm = verdicts(tree, base, key, [mid], flags)[0]
        if vm == vlo:
            lo = mid
        else:
            hi = mid
    return lo, hi


SCANS = [
    ("2_PEFT_gamma3_const_1", "EFTGamma30", np.linspace(-2, 2, 41), ["positivity", "minkowski", "mass", "mass_math"]),
    ("2_PEFT_Omega_const_1", "EFTOmega0", np.linspace(-2, 2, 41), ["positivity", "minkowski", "mass", "mass_math"]),
    ("2_PEFT_Omega_power_1", "EFTOmega0", np.linspace(-2, 2, 41), ["positivity", "minkowski", "mass", "mass_math"]),
    ("2_PEFT_gamma1_const_1", "EFTGamma10", np.linspace(-2, 2, 41), ["mass", "mass_math"]),
    ("4_Des_F_R_LCDM_1", "EFTB0", np.concatenate([-np.logspace(1, -6, 15), np.logspace(-6, 1, 15)]),
     ["positivity", "minkowski", "mass", "mass_math"]),
    ("5_Horava_eta_1", "Horava_eta", np.concatenate([-np.logspace(-1, -6, 11), np.logspace(-6, -1, 11)]),
     ["positivity", "minkowski", "mass", "mass_math"]),
    ("008p7-ext-galileon", "ExtendedGalileon_S", np.linspace(0.1, 5.0, 50), ["ini", "positivity", "minkowski", "mass",
                                                                             "mass_math"]),
    ("008p7-ext-galileon", "ExtendedGalileon_q", np.linspace(0.05, 2.0, 40), ["ini", "positivity", "minkowski", "mass",
                                                                              "mass_math"]),
]


def main():
    out_json, tb, tf = sys.argv[1:4]
    only = sys.argv[4:] or None
    report = []
    for model, key, grid, cfgs in SCANS:
        if only and model not in only:
            continue
        base = PAR / f"{model}.ini" if (PAR / f"{model}.ini").exists() else AUD / f"{model}.ini"
        for cfg in cfgs:
            flags = FLAGS[cfg]
            vb = verdicts(tb, base, key, grid, flags)
            vf = verdicts(tf, base, key, grid, flags)
            row = {"model": model, "key": key, "cfg": cfg, "grid": grid.tolist(), "buggy": "".join(vb),
                   "fixed": "".join(vf), "boundaries": []}
            print(f"{model:24s} {key:11s} {cfg:10s} buggy {''.join(vb)}\n{'':48s}fixed {''.join(vf)}", flush=True)
            for i in range(len(grid) - 1):
                if vb[i] != vb[i + 1] or vf[i] != vf[i + 1]:
                    bb = bisect(tb, base, key, grid[i], grid[i + 1], vb[i], flags) if vb[i] != vb[i + 1] else None
                    bf = bisect(tf, base, key, grid[i], grid[i + 1], vf[i], flags) if vf[i] != vf[i + 1] else None
                    row["boundaries"].append({"interval": [grid[i], grid[i + 1]], "buggy": bb, "fixed": bf})
                    print(f"      boundary in [{grid[i]:.4g},{grid[i + 1]:.4g}]: buggy {bb}  fixed {bf}", flush=True)
            report.append(row)
            Path(out_json).write_text(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()

"""Stability-verdict scan for finding pade-fourth-derivative: A0 vs one-line fix (work/verify/pade_fix).

Pure EFT (LCDM background) with Omega(a) = (N0 + N1 (a-a0) + N2 (a-a0)^2) / (1 + D1 (a-a0)) (Pade flag 6).
EFTOmegaPPPP enters EFTcdotdot / EFTLambdadotdot, i.e. only the mass-stability (EFT_mass_stability) and the
Minkowski-limit tachyon (EFT_minkowski_limit) conditions. Each tree x config runs in its own subprocess
(holding one harness CPU slot); the verdict is taken from camb.get_background(no_thermo=True).

usage: pade_stab_scan.py [npts_per_axis]
"""

import itertools
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

W = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify")
OUT = W / "pade_scan"
TREES = {"A0": str(harness.tree_path("A0")), "FIX": str(W / "pade_fix")}
INI = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"
CONFIGS = {
    "mass": {"EFT_mass_stability": "T", "EFT_minkowski_limit": "F", "EFT_positivity_bounds": "F"},
    "mink": {"EFT_mass_stability": "F", "EFT_minkowski_limit": "T", "EFT_positivity_bounds": "F"},
}


def point_overrides(p):
    a0, n0, n1, n2, d1 = p
    return {"PureEFTmodelOmega": "6", "EFTOmega_Pade_order_N": "2", "EFTOmega_Pade_order_D": "1",
            "EFTOmegaa0": repr(a0), "EFTOmegaN0": repr(n0), "EFTOmegaN1": repr(n1), "EFTOmegaN2": repr(n2),
            "EFTOmegaD1": repr(d1), "EFTCAMB_feedback_level": "0"}


def worker(cfg):
    sys.path.insert(0, cfg["tree"])
    import camb

    assert Path(camb.__file__).resolve().is_relative_to(Path(cfg["tree"]).resolve())
    wdir = Path(cfg["wdir"])
    res = []
    for i, p in enumerate(cfg["points"]):
        ov = dict(CONFIGS[cfg["config"]])
        ov.update(point_overrides(p))
        ini = wdir / f"p{i}.ini"
        ini.write_text(f"DEFAULT({INI})\n" + "".join(f"{k} = {v}\n" for k, v in ov.items()))
        try:
            cp = camb.read_ini(str(ini), no_validate=True)
            camb.get_background(cp, no_thermo=True)
            v = "ok"
        except Exception as e:  # noqa: BLE001
            v = str(e)[:80]
        res.append(v)
    Path(cfg["out"]).write_text(json.dumps(res))


def run(tree, config, points):
    wdir = OUT / f"{tree}_{config}"
    wdir.mkdir(parents=True, exist_ok=True)
    cfg = {"tree": TREES[tree], "config": config, "points": points, "wdir": str(wdir), "out": str(wdir / "v.json")}
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        r = subprocess.run([sys.executable, __file__, "_worker", json.dumps(cfg)], env=env, cwd=wdir,
                           capture_output=True, text=True, timeout=1500, preexec_fn=harness._nice)
    if r.returncode != 0:
        print(tree, config, "rc", r.returncode, r.stderr[-2000:])
        return None
    return json.loads((wdir / "v.json").read_text())


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    import numpy as np

    pts = []
    # family 1: a0 = 0, Omega(0) = 0
    for n1, n2, d1 in itertools.product(np.linspace(-0.2, 0.2, n), np.linspace(-0.2, 0.2, n),
                                        np.linspace(-0.6, 2.0, n)):
        pts.append((0.0, 0.0, float(n1), float(n2), float(d1)))
    # family 2: a0 = 1 (expansion around today), Omega(1) = N0
    for n0, n1, d1 in itertools.product(np.linspace(-0.1, 0.1, n), np.linspace(-0.3, 0.3, n),
                                        np.linspace(-0.6, 0.6, n)):
        pts.append((1.0, float(n0), float(n1), 0.0, float(d1)))
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(t, c) for c in CONFIGS for t in TREES]
    with ThreadPoolExecutor(4) as ex:
        out = dict(zip(jobs, ex.map(lambda j: run(j[0], j[1], pts), jobs)))
    summary = {"npoints": len(pts)}
    flips = []
    for c in CONFIGS:
        a, b = out[("A0", c)], out[("FIX", c)]
        if a is None or b is None:
            summary[c] = "worker failed"
            continue
        na, nb = sum(v == "ok" for v in a), sum(v == "ok" for v in b)
        diff = [i for i in range(len(pts)) if (a[i] == "ok") != (b[i] == "ok")]
        summary[c] = {"A0_ok": na, "FIX_ok": nb, "n_flipped": len(diff),
                      "A0_stable_FIX_unstable": sum(a[i] == "ok" for i in diff),
                      "A0_unstable_FIX_stable": sum(b[i] == "ok" for i in diff),
                      "errors_A0": sorted(set(a) - {"ok"}), "errors_FIX": sorted(set(b) - {"ok"})}
        flips += [{"config": c, "point(a0,N0,N1,N2,D1)": pts[i], "A0": a[i], "FIX": b[i]} for i in diff]
    summary["flips_first"] = flips[:20]
    print(json.dumps(summary, indent=1))
    (OUT / f"summary_n{n}.json").write_text(json.dumps({"summary": summary, "flips": flips, "points": pts}, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "_worker":
        worker(json.loads(sys.argv[2]))
    else:
        main()

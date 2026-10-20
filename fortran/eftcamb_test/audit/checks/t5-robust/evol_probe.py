"""Track 5: inspect pi(a) / pidot(a) around the EFT perturbation turn-on for one model.

Runs the model (default settings plus optional overrides) with get_scale_evolution at the requested k values and
prints, for a window of scale factors around the turn-on, pi, pidot and the quasi-static attractor
pi_QS = -h0 E / (C + k^2 D1 + k^4 D2) built from the cached pi-equation coefficients.

Usage: python evol_probe.py MODEL OUTNAME [--set key=value ...] [--k 0.01,0.1,1] [--tree A0] [--window 0.3,3]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t5-robust"


def probe(npz: Path, a_on: float, window=(0.5, 4.0), step=1, h0_mpc=70 / 299792.458):
    d = harness.load(npz)
    # the cache fields "a" and "k" overwrite the harness evol_a / evol_k arrays: both are (nk, na)
    a = d["evol_a"][0] if d["evol_a"].ndim == 2 else d["evol_a"]
    ks = d["evol_k"][:, 0] if d["evol_k"].ndim == 2 else d["evol_k"]
    h0 = h0_mpc
    lines = []
    for ik, k in enumerate(ks):
        pi, pid = d["evol_pi"][ik], d["evol_pidot"][ik]
        C, D1, D2, E = (d[f"evol_EFTpi{x}"][ik] for x in ("C", "D1", "D2", "E"))
        den = C + k**2 * D1 + k**4 * D2
        with np.errstate(divide="ignore", invalid="ignore"):
            piqs = np.where(den != 0, -h0 * E / den, np.nan)
        lines.append(f"k={k:g}")
        sel = np.where((a > a_on * window[0]) & (a < a_on * window[1]))[0][::step]
        for i in sel:
            lines.append(f"  a={a[i]:.4e} pi={pi[i]: .5e} pidot={pid[i]: .5e} piQS={piqs[i]: .5e} "
                         f"pi/piQS={pi[i] / piqs[i] if piqs[i] else np.nan: .5f}")
    return "\n".join(lines)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("model")
    p.add_argument("outname")
    p.add_argument("--set", nargs="*", default=[])
    p.add_argument("--k", default="0.01,0.1,1")
    p.add_argument("--tree", default="A0")
    p.add_argument("--window", default="0.5,4")
    p.add_argument("--aon", type=float, default=0.0)
    p.add_argument("--step", type=int, default=1)
    a = p.parse_args()
    ini = harness.PARAMS_DIR / f"{a.model}.ini"
    if not ini.exists():
        ini = Path(a.model)
    out = WORK / "evol" / a.outname / f"{ini.stem}.npz"
    ov = dict(s.split("=", 1) for s in a.set)
    if not out.exists():
        st = harness.run_model(a.tree, ini, out, overrides=ov, evol_k=[float(x) for x in a.k.split(",")])
        print(st["status"], st["elapsed"])
    keys = harness.parse_ini(ini)
    keys.update(ov)
    a_on = a.aon or float(keys.get("EFTCAMB_turn_on_time", "1e-2"))
    w = tuple(float(x) for x in a.window.split(","))
    h0 = float(keys.get("hubble", "70")) / 299792.458
    print(probe(out, a_on, w, a.step, h0))


if __name__ == "__main__":
    main()

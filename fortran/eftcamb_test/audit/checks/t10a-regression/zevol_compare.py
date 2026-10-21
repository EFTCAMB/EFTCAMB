"""Track 10a: compare CAMB's own perturbation evolution (CAMBdata.get_redshift_evolution) of A0 and B1.

For each model the transfer variables (Weyl potential, delta_cdm, T/E sources, ...) at fixed k are computed in both
trees with identical settings (IntTolBoost=2, optional omnuh2=0) and compared as max |A0-B1| / max|B1| over z in a
range.  A model is also compared to the GR run of the same tree to show the size of the EFT effect.

Usage: python zevol_compare.py OUTDIR MODEL [MODEL ...] [--nu0] [--ab 1]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

KS = [1e-3, 3e-3, 1e-2, 3e-2, 0.1]
Z = np.concatenate([np.linspace(0, 30, 121), np.logspace(np.log10(31), np.log10(1100), 120)])
VARS = ["delta_cdm", "delta_baryon", "delta_photon", "Weyl", "etak", "v_newtonian_cdm", "T_source", "E_source",
        "lens_potential_source"]


def worker(tree: str, ini: str, out: str, ov: dict) -> None:
    sys.path.insert(0, tree)
    import camb

    wrapper = Path(out).with_suffix(".wrapper.ini")
    wrapper.write_text(f"DEFAULT({ini})\n" + "".join(f"{k} = {v}\n" for k, v in ov.items())
                       + "EFTCAMB_write_background = F\n")
    cp = camb.read_ini(str(wrapper), no_validate=True)
    keys = harness.parse_ini(str(wrapper))
    cp.Accuracy.IntTolBoost = 2.0
    cp.SourceWindows = []  # number-count windows are irrelevant here
    if hasattr(cp, "EFTCAMB") and cp.EFTCAMB.EFTflag == 0 and int(keys.get("EFTflag", "0")) != 0:
        cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)
    r = camb.get_results(cp)
    ev = r.get_redshift_evolution(KS, Z, VARS)
    np.savez(out, ev=np.asarray(ev), k=np.asarray(KS), z=Z, omnuh2=cp.omnuh2)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("outdir")
    p.add_argument("models", nargs="*")
    p.add_argument("--nu0", action="store_true")
    p.add_argument("--worker", nargs=4)
    a = p.parse_args()
    if a.worker:
        tree, ini, out, ov = a.worker
        worker(tree, ini, out, json.loads(ov))
        return
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ov = {"omnuh2": "0.0"} if a.nu0 else {}
    models = ["1_EFT_GR"] + [m for m in a.models if m != "1_EFT_GR"]
    for m in models:
        for t in ("A0", "B1"):
            f = outdir / f"{t}_{m}.npz"
            if f.exists():
                continue
            env = dict(os.environ, OMP_NUM_THREADS="1")
            with harness.slot(1):
                subprocess.run(["nice", "-n", "10", sys.executable, __file__, str(outdir), "--worker",
                                str(Path(harness.AUDIT_ROOT) / "trees" / t), str(harness.PARAMS_DIR / f"{m}.ini"),
                                str(f), json.dumps(ov)], env=env, check=True, capture_output=True)
    d = {(t, m): np.load(outdir / f"{t}_{m}.npz") for m in models for t in ("A0", "B1")}
    late = Z <= 30
    for m in models:
        x, y = d[("A0", m)]["ev"], d[("B1", m)]["ev"]
        print(f"== {m} (omnuh2={float(d[('A0', m)]['omnuh2']):g}): max|A0-B1|/max|B1| for z<=30 | z>30")
        for ik, k in enumerate(KS):
            s = []
            for iv, v in enumerate(VARS):
                u, w = x[ik, :, iv], y[ik, :, iv]
                e1 = np.max(np.abs(u - w)[late]) / np.max(np.abs(w)[late])
                e2 = np.max(np.abs(u - w)[~late]) / np.max(np.abs(w)[~late])
                s.append(f"{v}={e1:.0e}|{e2:.0e}")
            print(f"  k={k:g}: " + " ".join(s))


if __name__ == "__main__":
    main()

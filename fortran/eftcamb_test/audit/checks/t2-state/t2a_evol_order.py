"""Track 2a: does get_scale_evolution of one k depend on the k evaluated before it in the same call (state kept in
the OpenMP-private EvolutionVars, reused across loop iterations on a thread)?  One process, OMP_NUM_THREADS=1:
evolution of k=K alone vs [K0, K] vs [K1, K]; prints the max relative difference per field.

Usage: python t2a_evol_order.py TREE MODEL_INI [K]"""

from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402


def main() -> None:
    tree, ini = str(harness.tree_path(sys.argv[1])), Path(sys.argv[2])
    kk = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
    os.environ["OMP_NUM_THREADS"] = "1"
    with harness.slot(1):
        os.nice(10)
        sys.path.insert(0, tree)
        import numpy as np

        import camb

        out = harness.AUDIT_ROOT / "work" / "t2-state" / "evolorder"
        w = harness._wrapper_ini(ini, out, ini.stem + "_order", {"get_transfer": "F", "get_tensor_cls": "F",
                                                                 "do_lensing": "F", "l_max_scalar": "200",
                                                                 "feedback_level": "0"})
        cp = camb.read_ini(str(w))
        r = camb.get_results(cp)
        a = np.logspace(-6, 0, 300)
        runs = {}
        for name, ks in (("alone", [kk]), ("after_1e-4", [1e-4, kk]), ("after_1", [1.0, kk]), ("alone_again", [kk])):
            _, arr = cp.EFTCAMB.get_scale_evolution(r, ks, a)
            runs[name] = {f: np.asarray(arr[f])[-1] for f in arr.dtype.names}
        skip = {"EFTcdotdot", "EFTGamma1PP"}  # never written (uninitialised, reported separately)
        for name in ("alone_again", "after_1e-4", "after_1"):
            worst = []
            for f, x in runs["alone"].items():
                if f in skip:
                    continue
                y = runs[name][f]
                if not np.array_equal(x, y, equal_nan=True):
                    sc = np.nanmax(np.abs(x))
                    worst.append((float(np.nanmax(np.abs(x - y)) / sc) if sc else float("inf"), f))
            worst.sort(reverse=True)
            print(f"{ini.stem} k={kk:g}: '{name}' vs 'alone': {len(worst)} fields differ; worst: "
                  + ", ".join(f"{f}={v:.1e}" for v, f in worst[:6]), flush=True)


if __name__ == "__main__":
    main()

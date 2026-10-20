"""Ionization history x_e(z) and visibility of a model in one tree (run in a subprocess per tree).

Usage: python xe_history.py TREE_PATH INI OUT_NPZ
"""

from __future__ import annotations

import json
import sys

import numpy as np


def main() -> None:
    tree, ini, out = sys.argv[1:4]
    sys.path.insert(0, tree)
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
    import camb
    import harness

    cp = camb.read_ini(ini, no_validate=True)
    keys = harness.parse_ini(ini)
    if hasattr(cp, "EFTCAMB") and cp.EFTCAMB.EFTflag == 0 and int(keys.get("EFTflag", "0")) != 0:
        cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)
    cp.WantCls = False
    cp.WantTransfer = False
    cp.SourceWindows = []  # upstream CAMB segfaults in get_background(no_thermo=False) with windows and no Cls
    r = camb.get_background(cp, no_thermo=False)
    z = np.concatenate([np.linspace(0, 40, 801), np.logspace(np.log10(41), np.log10(3000), 400)])
    ev = r.get_background_redshift_evolution(z, ["x_e", "visibility", "opacity"], format="array")
    d = r.get_derived_params()
    np.savez(out, z=z, x_e=ev[:, 0], vis=ev[:, 1], opac=ev[:, 2], meta=json.dumps({k: float(v) for k, v in d.items()}),
             zre=float(cp.Reion.redshift), tau=float(r.get_tau()) if hasattr(r, "get_tau") else -1.0)


if __name__ == "__main__":
    main()

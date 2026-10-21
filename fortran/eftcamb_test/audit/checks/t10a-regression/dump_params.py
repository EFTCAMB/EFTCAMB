"""Dump the CAMBparams produced by camb.read_ini for an ini in a given tree (run in a subprocess per tree).

Usage: python dump_params.py TREE_PATH INI [--after-results]
Prints a flattened key=value listing (Accuracy, Reion, Transfer, top-level scalars, etc.) as JSON.
"""

from __future__ import annotations

import ctypes
import json
import sys
from pathlib import Path


def flatten(obj, prefix="", depth=0, out=None):
    out = {} if out is None else out
    if depth > 3:
        return out
    fields = getattr(obj, "_fields_", None)
    names = []
    for cls in type(obj).__mro__:
        for f in getattr(cls, "_fields_", []) or []:
            names.append(f[0])
    for n in names:
        if n.startswith("_") or n in ("EFTCAMB", "eft_par_cache"):
            continue
        try:
            v = getattr(obj, n)
        except Exception as e:  # noqa: BLE001
            out[prefix + n] = f"<err {e!r}>"
            continue
        if isinstance(v, (int, float, bool, str)):
            out[prefix + n] = v
        elif isinstance(v, ctypes.Array):
            try:
                out[prefix + n] = [float(x) for x in v][:20]
            except Exception:  # noqa: BLE001
                out[prefix + n] = str(v)[:200]
        elif hasattr(v, "_fields_"):
            out[prefix + n + ".__class__"] = type(v).__name__
            flatten(v, prefix + n + ".", depth + 1, out)
        else:
            out[prefix + n] = str(v)[:200]
    return out


def main():
    tree, ini = sys.argv[1], sys.argv[2]
    sys.path.insert(0, tree)
    import camb

    cp = camb.read_ini(ini, no_validate=True)
    d = flatten(cp)
    if "--after-results" in sys.argv:
        r = camb.get_results(cp)
        d2 = flatten(r.Params, "RESULT.")
        d.update(d2)
    print(json.dumps(d, default=str))


if __name__ == "__main__":
    main()

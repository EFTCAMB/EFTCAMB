"""Default-kind (single precision) real literals in EFTCAMB code and EFTCAMB MOD blocks.

For every executable statement, find real literals without a kind suffix and without a 'd' exponent
(e.g. 0.1, 1.5e-3, 2.) and compute the relative change float32(v)/float64(v)-1.  Literals that are exactly
representable in single precision (0.5, 2., 0.25 ...) are harmless; the others are listed.
Usage: python check_real_literals.py [--all]
"""

from __future__ import annotations

import argparse
import collections
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import fparse  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))
LIT = re.compile(r"(?<![\w.])((?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?|\d+[eE][+-]?\d+)(?![\w.])(?!_)")


def strip_strings(s):
    return re.sub(r"'[^']*'|\"[^\"]*\"", "''", s)


def mod_ranges(path):
    rngs, start = [], None
    for i, ln in enumerate(open(path, errors="replace").read().split("\n"), 1):
        if "EFTCAMB MOD START" in ln:
            start = i
        elif ("EFTCAMB MOD END" in ln or "EFTCAMB MOD STOP" in ln) and start:
            rngs.append((start, i))
            start = None
    return rngs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    fs = sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/**/*.f90"), recursive=True))
    fs = [f for f in fs if "02_bundled" not in f]
    camb = [f for f in sorted(glob.glob(os.path.join(REPO, "fortran/*.f90"))) if mod_ranges(f)]
    per_file = collections.Counter()
    inexact_per_file = collections.Counter()
    rows = []
    for f in fs + camb:
        rng = mod_ranges(f) if f in camb else None
        for p in fparse.parse_file(f):
            for st in list(p.decls) + list(p.body):
                if rng is not None and not any(x <= st.line <= y for x, y in rng):
                    continue
                t = strip_strings(st.text)
                if re.match(r"^\s*(write|print|format)\b", t, re.I):
                    continue
                # strip format strings in write(*,'(...)')
                for m in LIT.finditer(t):
                    s = m.group(1)
                    v = float(s)
                    rel = 0.0 if v == 0 else abs(float(np.float32(v)) / v - 1)
                    rel_f = f"{rel:.1e}"
                    per_file[os.path.relpath(f, REPO)] += 1
                    if rel > 0:
                        inexact_per_file[os.path.relpath(f, REPO)] += 1
                        rows.append((os.path.relpath(f, REPO), st.line, p.name, s, rel_f, st.text[:150]))
                    elif a.all:
                        rows.append((os.path.relpath(f, REPO), st.line, p.name, s, rel_f, st.text[:150]))
    print("== per file: total kind-less real literals / inexact in single precision")
    for k in sorted(per_file):
        print(f"{per_file[k]:5d} {inexact_per_file[k]:5d}  {k}")
    print("TOTAL", sum(per_file.values()), "inexact", sum(inexact_per_file.values()))
    print("== inexact literals")
    for r in rows:
        print(f"{r[0]}:{r[1]}: {r[2]}: {r[3]} (rel {r[4]}) | {r[5]}")


if __name__ == "__main__":
    main()

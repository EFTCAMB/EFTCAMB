"""Copy-paste heuristic 2: assignments to DIFFERENT targets with an IDENTICAL (whitespace-normalised) right-hand side
of non-trivial length, within one file of the EFTCAMB sources or the EFTCAMB MOD blocks of fortran/*.f90.
Such pairs are either intentional aliases or a formula pasted under the wrong name (e.g. QS_sigma = <QS_mu formula>).

Usage: python check_duplicate_rhs.py [--min-len N]
"""

from __future__ import annotations

import collections
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import fparse  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))


def main():
    min_len = 60
    if "--min-len" in sys.argv:
        min_len = int(sys.argv[sys.argv.index("--min-len") + 1])
    fs = sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/**/*.f90"), recursive=True))
    fs = [f for f in fs if "02_bundled" not in f]
    n = 0
    for f in fs:
        by_rhs = collections.defaultdict(list)
        for p in fparse.parse_file(f):
            for st in p.body:
                a = fparse.assignment(st.text)
                if not a:
                    continue
                _, lhs, rhs = a
                key = re.sub(r"\s+", "", rhs).lower()
                if len(key) < min_len:
                    continue
                by_rhs[key].append((st.line, p.name, re.sub(r"\s+", "", lhs).lower()))
        for key, occ in by_rhs.items():
            lhss = {o[2] for o in occ}
            if len(lhss) > 1:
                n += 1
                rel = os.path.relpath(f, REPO)
                locs = ", ".join(f"{o[0]}:{o[1]}:{o[2]}" for o in occ)
                print(f"{rel}: same RHS ({len(key)} chars) assigned to {sorted(lhss)} at {locs}")
    print(f"TOTAL {n}")


if __name__ == "__main__":
    main()

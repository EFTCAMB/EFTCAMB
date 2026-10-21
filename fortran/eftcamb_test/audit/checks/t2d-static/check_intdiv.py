"""Integer division / integer power hygiene scan over EFTCAMB code and EFTCAMB MOD blocks.

Flags in executable statements:
  * a/b where both operands are integer literals or integer-typed names (declared in the procedure),
  * integer literal raised to a negative integer power (10**(-1) == 0),
  * real ** (int/int) exponents.
Only statements whose assignment target is not an integer variable are reported (likely physics expressions).
Usage: python check_intdiv.py
"""

from __future__ import annotations

import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import fparse  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))
INT_LIT = r"(?<![\w.])\d+(?![\w.])"
OPND = r"(?<![\w.%])([A-Za-z_]\w*|\d+)(?![\w.(%])"


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
    fs = sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/**/*.f90"), recursive=True))
    fs = [f for f in fs if "02_bundled" not in f]
    camb = [f for f in sorted(glob.glob(os.path.join(REPO, "fortran/*.f90"))) if mod_ranges(f)]
    hits = 0
    for f in fs + camb:
        rng = mod_ranges(f) if f in camb else None
        for p in fparse.parse_file(f):
            ints = set()
            for d in p.decls:
                tl = d.text.lower().strip()
                if tl.startswith("integer"):
                    _, ents = fparse.decl_entities(d.text)
                    ints |= {n for n, _ in ents}
            for st in p.body:
                if rng is not None and not any(a <= st.line <= b for a, b in rng):
                    continue
                t = strip_strings(st.text)
                tl = t.lower()
                if re.match(r"^\s*(write|print|format|call\s+mpistop)\b", tl):
                    continue
                asg = fparse.assignment(t)
                if asg and asg[0] in ints:
                    continue
                if re.match(r"^\s*do\s+\w+\s*=", tl):
                    continue
                reasons = []
                for m in re.finditer(OPND + r"\s*/\s*" + OPND, tl):
                    a, b = m.group(1), m.group(2)
                    ia = a.isdigit() or a in ints
                    ib = b.isdigit() or b in ints
                    # skip a/b inside array section or index context: crude check for preceding '(' with ':'
                    if ia and ib:
                        reasons.append(f"intdiv {a}/{b}")
                for m in re.finditer(r"(?<![\w.])(\d+)\s*\*\*\s*\(?\s*-\s*\d+", tl):
                    reasons.append(f"int**negint {m.group(0)}")
                for m in re.finditer(r"\*\*\s*\(\s*(\d+)\s*/\s*(\d+)\s*\)", tl):
                    reasons.append(f"pow(int/int) {m.group(0)}")
                if reasons:
                    hits += 1
                    print(f"{os.path.relpath(f, REPO)}:{st.line}: {p.name}: {'; '.join(reasons)} | {st.text[:140]}")
    print("TOTAL", hits)


if __name__ == "__main__":
    main()

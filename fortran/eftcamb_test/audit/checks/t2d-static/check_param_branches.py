"""Copy-paste check: ParameterNames / ParameterNamesLatex / ParameterValues of each model must walk the same
(bound, offset, sub-object) sequence.  Extracts, per procedure, the ordered list of
   (else) if (... i <= BOUND) ... if (i - OFFSET == j) call self%OBJ%<method>(j, ...)
triples and compares the three procedures of the same model.  Also checks the offset of each branch equals the
bound of the previous branch.
Usage: python check_param_branches.py
"""

from __future__ import annotations

import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import fparse  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))


def triples(p):
    out = []
    bound = None
    for st in p.body:
        tl = st.text.lower()
        m = re.match(r"^(else\s+)?if\s*\((.*)\)\s*then$", tl)
        if m:
            b = re.findall(r"\bi\s*<=\s*(\w+)", m.group(2))
            if b:
                bound = b[-1]
            continue
        m = re.search(r"if\s*\(\s*i\s*-\s*(\w+)\s*==\s*j\s*\)\s*call\s+self%(\w+)%(\w+)", tl)
        if m:
            out.append((bound, m.group(1), m.group(2)))
            continue
        m = re.search(r"if\s*\(\s*i\s*==\s*j\s*\)\s*call\s+self%(\w+)%(\w+)", tl)
        if m:
            out.append((bound, "0", m.group(1)))
    return out


def main():
    fs = sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/0[78]f_*/*.f90")))
    issues = 0
    for f in fs:
        procs = {p.name.lower(): p for p in fparse.parse_file(f)}
        groups = {}
        for n, p in procs.items():
            for suf, key in (("parameternameslatex", "latex"), ("parameternames", "names"), ("parametervalues", "values")):
                if n.endswith(suf):
                    groups[key] = triples(p)
                    break
        if not groups:
            continue
        rel = os.path.relpath(f, REPO)
        ref = groups.get("names")
        for key, tr in groups.items():
            if ref is not None and tr != ref:
                issues += 1
                print(f"{rel}: {key} differs from names:")
                for x, y in zip(ref, tr):
                    if x != y:
                        print(f"     names {x}  vs  {key} {y}")
                if len(ref) != len(tr):
                    print(f"     lengths {len(ref)} vs {len(tr)}")
            prev = "0"
            for bnd, off, obj in tr:
                if off != prev and off != "0":
                    issues += 1
                    print(f"{rel}: {key}: branch {obj} (bound {bnd}) uses offset {off}, previous bound {prev}")
                prev = bnd
        print(f"{rel}: branches names={len(groups.get('names', []))} latex={len(groups.get('latex', []))} "
              f"values={len(groups.get('values', []))}")
    print("TOTAL issues", issues)


if __name__ == "__main__":
    main()

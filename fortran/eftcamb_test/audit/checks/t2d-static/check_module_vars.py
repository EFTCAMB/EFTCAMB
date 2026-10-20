"""List module-level (non-parameter) variables in EFTCAMB sources and where they are written."""

from __future__ import annotations

import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import fparse  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))


def module_level_decls(path):
    stmts = fparse.logical_statements(path)
    out = []
    depth_proc = 0
    in_type = 0
    in_iface = 0
    in_module = False
    for st in stmts:
        tl = st.text.lower().strip()
        if re.match(r"^module\s+(?!procedure)\w+\s*$", tl) or re.match(r"^submodule\b", tl):
            in_module = True
            continue
        if re.match(r"^end\s+(sub)?module\b", tl):
            in_module = False
            continue
        if re.match(r"^(abstract\s+)?interface\b", tl):
            in_iface += 1
            continue
        if re.match(r"^end\s+interface\b", tl):
            in_iface -= 1
            continue
        if in_iface:
            continue
        if re.match(r"^type\s*(,[^:]*)?(::)?\s*\w+\s*$", tl) and not tl.startswith("type("):
            in_type += 1
            continue
        if re.match(r"^end\s*type\b", tl):
            in_type -= 1
            continue
        if in_type:
            continue
        if tl == "contains" and depth_proc == 0:
            in_module = False  # past spec part
            continue
        if in_module and depth_proc == 0 and fparse.TYPE_RE.match(tl) and "function" not in tl.split("::")[0]:
            spec, ents = fparse.decl_entities(st.text)
            if "parameter" in spec.lower():
                continue
            for n, v in ents:
                out.append((st.line, n, spec.strip(), v))
    return out


def main():
    fs = sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/**/*.f90"), recursive=True))
    allv = {}
    for f in fs:
        for line, n, spec, v in module_level_decls(f):
            allv[n] = (os.path.relpath(f, REPO), line, spec, v)
            print(f"{os.path.relpath(f, REPO)}:{line}: {spec} :: {n} {'= ' + v if v else ''}")
    print("TOTAL", len(allv))
    # writes
    print("\n--- writes to module variables inside procedures ---")
    for f in fs + sorted(glob.glob(os.path.join(REPO, "fortran/*.f90"))):
        for p in fparse.parse_file(f):
            for st in p.body:
                a = fparse.assignment(st.text)
                if a and a[0] in allv:
                    print(f"{os.path.relpath(f, REPO)}:{st.line}: {p.name}: {a[1]} = ...  (decl {allv[a[0]][0]}:{allv[a[0]][1]})")


if __name__ == "__main__":
    main()

"""List locals initialised in their declaration (implicit SAVE) and explicit SAVE inside procedures.

Usage: python check_implicit_save.py [--json out.json]
Scans fortran/eftcamb/**/*.f90 and the CAMB files carrying EFTCAMB MOD blocks.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import fparse  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))


def files() -> list[str]:
    fs = sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/**/*.f90"), recursive=True))
    for f in sorted(glob.glob(os.path.join(REPO, "fortran/*.f90"))):
        if "EFTCAMB MOD START" in open(f, errors="replace").read():
            fs.append(f)
    return fs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json")
    a = ap.parse_args()
    rows = []
    for f in files():
        for p in fparse.parse_file(f):
            if p.in_interface:
                continue
            for d in p.decls:
                tl = d.text.lower()
                if re.match(r"^\s*save\b", tl):
                    rows.append(dict(file=os.path.relpath(f, REPO), line=d.line, proc=p.name, var="<SAVE stmt>", init=d.text))
                    continue
                if not ffparse_is_typedecl(tl):
                    continue
                spec, ents = fparse.decl_entities(d.text)
                if "parameter" in spec.lower():
                    continue
                for n, v in ents:
                    if v is not None:
                        rows.append(
                            dict(
                                file=os.path.relpath(f, REPO),
                                line=d.line,
                                proc=p.name,
                                module=p.module,
                                var=n,
                                init=v,
                                spec=spec.strip(),
                                save_attr="save" in spec.lower(),
                            )
                        )
    for r in rows:
        print(f"{r['file']}:{r['line']}: {r['proc']}: {r['var']} = {r['init']}")
    print(f"TOTAL {len(rows)}")
    if a.json:
        json.dump(rows, open(a.json, "w"), indent=1)


def ffparse_is_typedecl(tl: str) -> bool:
    return bool(fparse.TYPE_RE.match(tl))


if __name__ == "__main__":
    main()

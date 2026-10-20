"""Scan EFTCAMB Fortran for subroutines with a `logical, intent(out) :: success` dummy that never assign .true. to it.

Heuristic (text-level): per subroutine body (excluding contained procedures), list dummies named success declared
intent(out) and report those with no `success = .true.` assignment (they are only defined on an error path).
Usage: python scan_intent_out_success.py REPO
"""

import re
import sys
from pathlib import Path

repo = Path(sys.argv[1])
pat_start = re.compile(r"^\s*(?:recursive\s+|pure\s+|elemental\s+)*subroutine\s+(\w+)", re.I)
pat_end = re.compile(r"^\s*end\s+subroutine", re.I)
for f in sorted((repo / "fortran" / "eftcamb").rglob("*.f90")):
    lines = f.read_text(errors="replace").splitlines()
    stack: list[tuple[str, int, list[str]]] = []
    for i, ln in enumerate(lines, 1):
        code = ln.split("!", 1)[0]
        m = pat_start.match(code)
        if m:
            stack.append((m.group(1), i, []))
            continue
        if pat_end.match(code) and stack:
            name, start, body = stack.pop()
            txt = "\n".join(body)
            if re.search(r"logical[^:\n]*intent\s*\(\s*out\s*\)[^\n]*::\s*success\b", txt, re.I) and not re.search(
                r"^\s*success\s*=\s*\.true\.", txt, re.I | re.M
            ):
                sets = [b.strip() for b in body if re.match(r"^\s*success\s*=", b, re.I)]
                print(f"{f.relative_to(repo)}:{start} {name}: no 'success = .true.'; assignments: {sets[:3]}")
            continue
        if stack:
            stack[-1][2].append(code)

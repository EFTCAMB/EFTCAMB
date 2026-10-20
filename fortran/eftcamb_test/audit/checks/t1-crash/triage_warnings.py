"""Parse a gfortran build log (built with -Wall -Wextra -Wconversion ...) and triage the warnings.

Usage: python triage_warnings.py BUILD_LOG [--scope eft|mods|all] [--cat CATEGORY] [--list]

Scope 'eft' = files under fortran/eftcamb/, 'mods' = CAMB files that contain EFTCAMB MOD blocks (warnings inside the
blocks only, using line ranges from the source file), 'all' = everything.
"""

from __future__ import annotations

import argparse
import collections
import json
import re
from pathlib import Path

LOC = re.compile(r"^(/[^:]+\.(?:f90|f|F90|c)):(\d+):(\d+):\s*$")
LOC2 = re.compile(r"^(/[^:]+\.(?:f90|f|F90|c)):(\d+):(\d+)(?:-\d+)?:\s*(Warning|Error):\s*(.*)$")
MSG = re.compile(r"^(Warning|Error):\s*(.*?)(?:\s*\[(-W[\w=-]+)\])?$")


def parse(log: Path) -> list[dict]:
    lines = log.read_text(errors="replace").splitlines()
    out, cur = [], None
    for ln in lines:
        m = LOC.match(ln)
        if m:
            cur = (m.group(1), int(m.group(2)), int(m.group(3)))
            continue
        m2 = LOC2.match(ln)
        if m2:
            f, l, c, kind, msg = m2.groups()
            mm = MSG.match(f"{kind}: {msg}")
            out.append({"file": f, "line": int(l), "col": int(c), "kind": kind, "msg": mm.group(2),
                        "cat": mm.group(3) or "none"})
            continue
        m = MSG.match(ln)
        if m and cur:
            out.append({"file": cur[0], "line": cur[1], "col": cur[2], "kind": m.group(1), "msg": m.group(2),
                        "cat": m.group(3) or "none"})
            cur = None
    # dedupe (same location+message)
    seen, res = set(), []
    for w in out:
        key = (w["file"], w["line"], w["col"], w["msg"])
        if key not in seen:
            seen.add(key)
            res.append(w)
    return res


def mod_ranges(path: Path) -> list[tuple[int, int]]:
    rng, start = [], None
    for i, ln in enumerate(path.read_text(errors="replace").splitlines(), 1):
        if "EFTCAMB MOD START" in ln or "EFTCAMB MOD start" in ln:
            start = i
        elif ("EFTCAMB MOD END" in ln or "EFTCAMB MOD end" in ln) and start:
            rng.append((start, i))
            start = None
    return rng


def in_scope(w: dict, scope: str, cache: dict) -> bool:
    f = w["file"]
    if scope == "all":
        return True
    if "/fortran/eftcamb/" in f:
        return scope in ("eft", "eftmods")
    if scope in ("mods", "eftmods") and "/fortran/" in f and "/forutils/" not in f:
        if f not in cache:
            cache[f] = mod_ranges(Path(f))
        return any(a <= w["line"] <= b for a, b in cache[f])
    return False


def short(f: str) -> str:
    i = f.find("/fortran/")
    return f[i + 1 :] if i >= 0 else f


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("log")
    p.add_argument("--scope", default="eftmods")
    p.add_argument("--cat", default="")
    p.add_argument("--grep", default="")
    p.add_argument("--list", action="store_true")
    p.add_argument("--json", default="")
    a = p.parse_args()
    ws = parse(Path(a.log))
    cache: dict = {}
    ws = [w for w in ws if in_scope(w, a.scope, cache)]
    if a.cat:
        ws = [w for w in ws if w["cat"] == a.cat]
    if a.grep:
        ws = [w for w in ws if re.search(a.grep, w["msg"])]
    cats = collections.Counter(w["cat"] for w in ws)
    print(f"total={len(ws)}")
    for c, n in cats.most_common():
        print(f"  {c:40s} {n}")
    byfile = collections.Counter(short(w["file"]) for w in ws)
    if not a.list:
        print("by file:")
        for f, n in byfile.most_common(60):
            print(f"  {f:70s} {n}")
    else:
        for w in sorted(ws, key=lambda w: (w["file"], w["line"])):
            print(f"{short(w['file'])}:{w['line']}:{w['col']}: [{w['cat']}] {w['msg']}")
    if a.json:
        Path(a.json).write_text(json.dumps(ws, indent=1))


if __name__ == "__main__":
    main()

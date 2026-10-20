"""Parse a gfortran warnings build log into records; filter to EFTCAMB code (fortran/eftcamb and EFTCAMB MOD blocks).

Usage: python parse_warnings.py LOG TREE [--kind KIND] [--json OUT] [--summary]
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import re

HDR = re.compile(r"^(?P<file>[^\s:][^:]*\.f90):(?P<line>\d+):(?P<col>\d+):\s*$")
WARN = re.compile(r"^(Warning|Error): (?P<msg>.*?)(\s*\[(?P<kind>-W[\w-]+)\])?\s*$")


def mod_ranges(path):
    rngs = []
    start = None
    try:
        lines = open(path, errors="replace").read().split("\n")
    except OSError:
        return rngs
    for i, ln in enumerate(lines, 1):
        if "EFTCAMB MOD START" in ln:
            start = i
        elif "EFTCAMB MOD END" in ln or "EFTCAMB MOD STOP" in ln:
            if start:
                rngs.append((start, i))
            start = None
    return rngs


def parse(log, tree):
    recs = []
    cur = None
    src_lines = []
    for ln in open(log, errors="replace"):
        ln = ln.rstrip("\n")
        m = HDR.match(ln)
        if m:
            cur = dict(file=m.group("file"), line=int(m.group("line")), col=int(m.group("col")))
            src_lines = []
            continue
        if cur is not None:
            w = WARN.match(ln)
            if w:
                cur["msg"] = w.group("msg")
                cur["kind"] = w.group("kind") or ""
                cur["src"] = "\n".join(src_lines)
                recs.append(cur)
                cur = None
            else:
                src_lines.append(ln)
    # normalise paths
    seen = set()
    out = []
    mr_cache = {}
    for r in recs:
        f = r["file"]
        if not os.path.isabs(f):
            f = os.path.normpath(os.path.join(tree, "fortran", f))
        rel = os.path.relpath(f, tree)
        r["rel"] = rel
        key = (rel, r["line"], r["col"], r["kind"], r["msg"])
        if key in seen:
            continue
        seen.add(key)
        if rel.startswith("fortran/eftcamb/"):
            r["scope"] = "bundled" if "02_bundled" in rel else "eftcamb"
        else:
            if rel not in mr_cache:
                mr_cache[rel] = mod_ranges(f)
            r["scope"] = "mod" if any(a <= r["line"] <= b for a, b in mr_cache[rel]) else "camb"
        out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("tree")
    ap.add_argument("--json")
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()
    recs = parse(a.log, a.tree)
    if a.json:
        json.dump(recs, open(a.json, "w"), indent=1)
    if a.summary:
        c = collections.Counter((r["scope"], r["kind"]) for r in recs)
        for k, v in sorted(c.items()):
            print(f"{k[0]:8s} {k[1]:28s} {v}")


if __name__ == "__main__":
    main()

"""intent(out) dummies and function results that are never (or not on all simple paths) assigned.

Heuristic: a name counts as set when it appears as an assignment target, as an actual argument of a call
(may be set by the callee), in a read/allocate statement, or as the target of a pointer assignment.
Additionally reports a 'return' statement reached before the first set (straight-line order), which flags
early-return paths leaving the dummy/result undefined.
Usage: python check_intent_out.py
"""

from __future__ import annotations

import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import fparse  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))


def mod_ranges(path):
    rngs, start = [], None
    for i, ln in enumerate(open(path, errors="replace").read().split("\n"), 1):
        if "EFTCAMB MOD START" in ln:
            start = i
        elif ("EFTCAMB MOD END" in ln or "EFTCAMB MOD STOP" in ln) and start:
            rngs.append((start, i))
            start = None
    return rngs


def sets_name(st_text, name):
    tl = st_text.lower()
    asg = fparse.assignment(st_text)
    if asg and asg[0] == name:
        return True
    if re.match(r"^\s*(if\s*\(.*\)\s*)?call\b", tl) and re.search(r"[(,]\s*" + name + r"\b", tl):
        return True
    if re.match(r"^\s*(read|allocate|inquire)\b", tl) and re.search(r"\b" + name + r"\b", tl):
        return True
    if re.search(r"^\s*" + name + r"\s*(\(.*\))?\s*=>", tl):
        return True
    # function result used as argument to a subroutine call is covered above; also 'name%x = ' handled by assignment
    return False


def main():
    fs = sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/**/*.f90"), recursive=True))
    fs = [f for f in fs if "02_bundled" not in f]
    n_never = n_early = 0
    for f in fs:
        rel = os.path.relpath(f, REPO)
        for p in fparse.parse_file(f):
            if p.in_interface:
                continue
            outs = []
            for d in p.decls:
                tl = d.text.lower()
                if re.search(r"intent\s*\(\s*out\s*\)", tl):
                    _, ents = fparse.decl_entities(d.text)
                    outs += [e for e, _ in ents if e in p.args]
            if p.kind == "function" and p.result:
                outs.append(p.result)
            for name in outs:
                first_set = None
                first_ret = None
                for st in p.body:
                    tl = st.text.lower().strip()
                    if first_set is None and sets_name(st.text, name):
                        first_set = st.line
                    if first_ret is None and re.match(r"^(if\s*\(.*\)\s*)?return\b", tl):
                        first_ret = st.line
                if first_set is None:
                    n_never += 1
                    print(f"NEVER {rel}:{p.line}: {p.name}: '{name}' never assigned")
                elif first_ret is not None and first_ret < first_set:
                    n_early += 1
                    print(f"EARLY {rel}:{first_ret}: {p.name}: return before first set of '{name}' (line {first_set})")
    print("TOTAL never", n_never, "early-return", n_early)


if __name__ == "__main__":
    main()

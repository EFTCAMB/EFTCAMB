"""Extract EFTCAMB MOD START/END blocks from equations.f90 at HEAD and at eftcamb-pre-rebase,
pair them in order (with fuzzy matching on the START comment) and report whitespace-insensitive
differences of the block bodies.

Usage: python compare_blocks.py [outdir]
"""

import difflib
import json
import os
import re
import subprocess
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", ".."))
OUT = (
    sys.argv[1]
    if len(sys.argv) > 1
    else "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
    "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/t8a-equations"
)

START = re.compile(r"!\s*EFTCAMB MOD START", re.I)
END = re.compile(r"!\s*EFTCAMB MOD (END|STOP)", re.I)


def git_show(ref, path="fortran/equations.f90"):
    return subprocess.check_output(["git", "-C", REPO, "show", f"{ref}:{path}"], text=True).splitlines()


def blocks(lines):
    out = []
    cur = None
    sub = None
    for i, line in enumerate(lines, 1):
        if START.search(line):
            if cur is not None:
                out.append(dict(start=cur, end=None, head=lines[cur - 1].strip(), body=lines[cur:i - 1], sub=sub))
            cur = i
            sub = current_sub(lines, i)
        elif END.search(line) and cur is not None:
            out.append(dict(start=cur, end=i, head=lines[cur - 1].strip(), body=lines[cur:i - 1], sub=sub))
            cur = None
    return out


SUBRE = re.compile(r"^\s*(recursive\s+|module\s+)?(subroutine|function)\s+(\w+)", re.I)


def current_sub(lines, i):
    for j in range(i - 1, -1, -1):
        m = SUBRE.match(lines[j])
        if m and not lines[j].strip().lower().startswith("end"):
            return m.group(3)
    return "?"


def norm(body):
    res = []
    for line in body:
        s = line.strip()
        if not s:
            continue
        s = re.sub(r"\s+", "", s).lower()
        res.append(s)
    return res


def main():
    head = git_show("HEAD")
    pre = git_show("eftcamb-pre-rebase")
    bh = blocks(head)
    bp = blocks(pre)
    # pair by order within the same subroutine name, using sequence matching of (sub, normalized body)
    keys_h = [(b["sub"].lower(), "\n".join(norm(b["body"]))) for b in bh]
    keys_p = [(b["sub"].lower(), "\n".join(norm(b["body"]))) for b in bp]
    sm = difflib.SequenceMatcher(a=[k[0] for k in keys_p], b=[k[0] for k in keys_h], autojunk=False)
    pairs = []
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal" or (tag == "replace" and i2 - i1 == j2 - j1):
            for di in range(i2 - i1):
                pairs.append((i1 + di, j1 + di))
        else:
            for i in range(i1, i2):
                pairs.append((i, None))
            for j in range(j1, j2):
                pairs.append((None, j))
    report = []
    for ip, ih in pairs:
        rec = {}
        if ih is not None:
            b = bh[ih]
            rec.update(head_start=b["start"], head_end=b["end"], sub=b["sub"], head=b["head"])
        if ip is not None:
            b = bp[ip]
            rec.update(pre_start=b["start"], pre_end=b["end"], pre_sub=b["sub"], pre_head=b["head"])
        if ih is not None and ip is not None:
            a = norm(bp[ip]["body"])
            c = norm(bh[ih]["body"])
            rec["identical"] = a == c
            if a != c:
                rec["diff"] = list(difflib.unified_diff(a, c, lineterm="", n=0))
        report.append(rec)
    with open(os.path.join(OUT, "block_pairs.json"), "w") as f:
        json.dump(report, f, indent=1)
    print(f"HEAD blocks: {len(bh)}  pre-rebase blocks: {len(bp)}")
    for r in report:
        tag = "SAME" if r.get("identical") else ("DIFF" if "diff" in r else ("ONLY-HEAD" if "head" in r else "ONLY-PRE"))
        print(
            f"{tag:9s} HEAD {r.get('head_start', '-')!s:>5}-{r.get('head_end', '-')!s:<5} "
            f"PRE {r.get('pre_start', '-')!s:>5}-{r.get('pre_end', '-')!s:<5} "
            f"{r.get('sub', r.get('pre_sub'))}: {r.get('head', r.get('pre_head'))[:70]}"
        )
        if "diff" in r:
            for d in r["diff"][2:]:
                print("            " + d[:160])


if __name__ == "__main__":
    main()

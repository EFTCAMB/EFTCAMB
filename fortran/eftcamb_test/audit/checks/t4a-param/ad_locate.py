"""Track t4a-param: locate the FPE trap / SIGSEGV of every failing AD driver case with lldb (frame 0 and 1 + the
x being evaluated).  Reads the case files written by run_param_checks.py --build AD.

usage: python ad_locate.py [names...]
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
W = AUDIT / "work" / "t4a-param"
EXE = W / "driver_AD"
CASES = W / "AD" / "cases"


def locate(name: str) -> str:
    case = CASES / f"{name}.case"
    cmd = ["lldb", "--batch", "-o", "run", "-k", "bt 4", "-k", "frame select 1", "-k", "frame variable i",
           "-k", "quit", "--", str(EXE), str(case)]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        out = p.stdout + p.stderr
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        subprocess.run(["pkill", "-f", f"driver_AD {case}"])
    stop = re.findall(r"stop reason = ([^\n]*)", out)
    frames = re.findall(r"frame #(\d): 0x[0-9a-f]+ ([^`]+)`([^ ]+)(?: at ([^\n]+))?", out)
    ivar = re.findall(r"\(int\) i = (-?\d+)", out)
    xi = None
    if ivar:
        i = int(ivar[-1])
        lines = case.read_text().splitlines()
        xi = lines[7 + i - 1] if 0 < i <= len(lines) - 7 else None
    fr = "; ".join(f"#{n} {fn} {loc}" for n, _, fn, loc in frames[:3])
    return f"{name}: {stop[:1]} | {fr} | i={ivar[-1:] } x={xi}"


def main():
    names = sys.argv[1:]
    for n in names:
        print(locate(n), flush=True)


if __name__ == "__main__":
    main()

"""Symbolise raw gfortran backtraces in trap_run.py logs with atos (macOS).

Each log written by trap_run.py starts with "IMAGE_BASE <path> <hex>" lines; the raw "#n 0x..." frames after
"Program received signal" are resolved against <path>.dSYM (created with dsymutil if missing).

Usage: python symbolize.py LOG [LOG ...] [--frames 8] [--brief]
  --brief prints one line per log: model, first in-camblib frame file:line and function.
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path


def dsym_for(lib: str) -> str:
    d = Path(lib + ".dSYM")
    if not d.exists():
        subprocess.run(["dsymutil", lib, "-o", str(d)], capture_output=True)
    dw = d / "Contents" / "Resources" / "DWARF" / Path(lib).name
    return str(dw) if dw.exists() else lib


def symbolize(log: Path, frames: int) -> list[str]:
    txt = log.read_text(errors="replace")
    bases = dict((m.group(1), int(m.group(2), 16)) for m in re.finditer(r"IMAGE_BASE (\S+) (0x[0-9a-f]+)", txt))
    camb = next((p for p in bases if p.endswith("camblib.so")), None)
    i = txt.find("Backtrace for this error")
    if i < 0 or camb is None:
        return []
    addrs = re.findall(r"^#\d+\s+(0x[0-9a-f]+)", txt[i:], re.M)
    if not addrs:
        return []
    base = bases[camb]
    # keep frames that fall inside camblib (heuristic: within 64 MB above its base)
    inlib = [a for a in addrs if 0 <= int(a, 16) - base < 64 << 20][:frames]
    if not inlib:
        return []
    out = subprocess.run(["atos", "-o", dsym_for(camb), "-l", hex(base), *inlib], capture_output=True, text=True)
    return [ln.strip() for ln in out.stdout.splitlines() if ln.strip()]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("logs", nargs="+", type=Path)
    p.add_argument("--frames", type=int, default=8)
    p.add_argument("--brief", action="store_true")
    a = p.parse_args()
    for log in a.logs:
        fr = symbolize(log, a.frames)
        if a.brief:
            first = fr[0] if fr else "(no backtrace)"
            m = re.match(r"(\S+) \(in camblib.so\) \((.*)\)", first)
            print(f"{log.stem:40s} {m.group(2) if m else '':38s} {m.group(1) if m else first}")
        else:
            print(f"== {log}")
            for f in fr:
                print("   ", f)


if __name__ == "__main__":
    main()

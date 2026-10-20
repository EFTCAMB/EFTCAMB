"""Run a command while holding N machine-wide CPU slots (harness.slot), at nice 10, with OMP_NUM_THREADS=N.

Usage: python slotrun.py N [--timeout S] [--out FILE] -- cmd arg ...
Used to run the standalone repro scripts (which spawn their own CAMB child processes) under the audit's resource rules.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402


def main() -> int:
    args = sys.argv[1:]
    n = int(args.pop(0))
    timeout, out = 560.0, None
    while args and args[0] != "--":
        opt = args.pop(0)
        if opt == "--timeout":
            timeout = float(args.pop(0))
        elif opt == "--out":
            out = args.pop(0)
    cmd = args[1:]
    env = dict(os.environ, OMP_NUM_THREADS=str(n))
    with harness.slot(n):
        try:
            p = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=timeout,
                               preexec_fn=lambda: os.nice(10))
            text, rc = p.stdout + p.stderr, p.returncode
        except subprocess.TimeoutExpired as e:
            text, rc = (e.stdout or b"").decode(errors="replace") + "\nTIMEOUT", -999
    if out:
        Path(out).write_text(text)
    print(text)
    print(f"[slotrun rc={rc}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())

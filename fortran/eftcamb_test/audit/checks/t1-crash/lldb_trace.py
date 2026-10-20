"""Re-run one harness model under lldb with FP traps enabled and print the exact faulting location + backtrace.

The libgfortran backtrace printed from inside Python is imprecise for FP traps (it reports the caller frame);
lldb stops on the faulting instruction (EXC_BAD_INSTRUCTION on Apple arm64).

Usage: python lldb_trace.py TREE MODEL_INI OUTDIR [--set key=value ...] [--notrap]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

DRIVER = """
import json, sys
sys.path.insert(0, {here!r})
import trap_run
cfg = json.load(open(sys.argv[1]))
if not {notrap!r}:
    trap_run._install_traps(cfg["tree"])
trap_run.harness._worker(cfg)
"""


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("tree")
    p.add_argument("ini")
    p.add_argument("outdir")
    p.add_argument("--set", nargs="*", default=[])
    p.add_argument("--notrap", action="store_true")
    p.add_argument("--frames", type=int, default=14)
    a = p.parse_args()
    out = Path(a.outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    ini = Path(a.ini)
    over = dict(s.split("=", 1) for s in a.set)
    wrapper = harness._wrapper_ini(ini, out, ini.stem + "_lldb", over)
    cfg = {"tree": str(harness.tree_path(a.tree)), "ini": str(wrapper), "model_ini": str(ini.resolve()),
           "out": str(out / f"{ini.stem}_lldb.npz"), "evol_k": [], "python_set": {}}
    (out / "lldb_cfg.json").write_text(json.dumps(cfg))
    (out / "lldb_drv.py").write_text(DRIVER.format(here=str(HERE), notrap=a.notrap))
    cmd = ["lldb", "--batch", "-o", "process handle SIGILL SIGFPE SIGSEGV SIGBUS SIGABRT -s true -p false",
           "-o", "breakpoint set -n _gfortran_runtime_error_at", "-o", "breakpoint set -n _gfortran_runtime_error",
           "-o", "breakpoint set -n _gfortran_stop_string", "-o", "breakpoint set -n _gfortran_stop_numeric",
           "-o", "run", "-o", f"bt {a.frames}",
           "-k", f"bt {a.frames}", "-k", "frame variable", "--", sys.executable, str(out / "lldb_drv.py"),
           str(out / "lldb_cfg.json")]
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        r = subprocess.run(cmd, cwd=out, env=env, capture_output=True, text=True, preexec_fn=harness._nice,
                           timeout=3000)
    txt = r.stdout + r.stderr
    (out / f"{ini.stem}_lldb.txt").write_text(txt)
    keep, on = [], False
    for ln in txt.splitlines():
        if "stop reason" in ln or ln.startswith("(lldb) bt") or "frame #" in ln:
            on = True
        if on:
            keep.append(ln)
    print("\n".join(keep[:200]) if keep else txt[-3000:])


if __name__ == "__main__":
    main()

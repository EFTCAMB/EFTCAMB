"""Shared helpers for the track-9 (Python wrapper) audit checks.

Every snippet runs in its own subprocess that imports camb only from the chosen frozen tree, so that crashes
(error stop, segfaults) are isolated and reported as a status instead of killing the driver.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REPO = HERE.parents[4]
PARAMS = REPO / "fortran" / "eftcamb_test" / "parameters"
AUDIT = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
    "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
)
TREES = AUDIT / "trees"
WORK = AUDIT / "work" / "t9-wrapper"
PY = "/Users/marco/Software/miniforge3/bin/python"

PRELUDE = """
import sys, json, os
sys.path.insert(0, {tree!r})
import numpy as np
import camb
from pathlib import Path
assert Path(camb.__file__).resolve().is_relative_to(Path({tree!r}).resolve()), camb.__file__
PARAMS = Path({params!r})
WORK = Path({work!r})
def emit(**kw):
    print("@@RESULT@@" + json.dumps(kw, default=lambda o: repr(o)), flush=True)
"""


def tree(name: str) -> str:
    p = TREES / name
    return str(p if p.exists() else Path(name))


def run_snippet(code: str, tree_name: str = "A0", threads: int = 2, timeout: float = 900, cwd: str | None = None,
                env_extra: dict | None = None) -> dict:
    """Run python code in a subprocess importing camb from the tree. Returns status, emitted results and log tail."""
    src = PRELUDE.format(tree=tree(tree_name), params=str(PARAMS), work=str(WORK)) + textwrap.dedent(code)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), PYTHONHASHSEED="0")
    env.pop("PYTHONPATH", None)
    env.update(env_extra or {})
    cwd = cwd or str(WORK)
    Path(cwd).mkdir(parents=True, exist_ok=True)
    try:
        # machine-wide CPU slots + nice 10 (shared machine resource rules)
        with harness.slot(threads):
            proc = subprocess.run([PY, "-c", src], cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout,
                                  preexec_fn=lambda: os.nice(10))
        rc = proc.returncode
        out, err = proc.stdout, proc.stderr
        status = "ok" if rc == 0 else (f"signal:{signal.Signals(-rc).name}" if rc < 0 else f"error:rc={rc}")
    except subprocess.TimeoutExpired as e:
        rc, status = None, "timeout"
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        err = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or "")
    results = []
    for line in out.splitlines():
        if line.startswith("@@RESULT@@"):
            results.append(json.loads(line[len("@@RESULT@@"):]))
    merged = {}
    for r in results:
        merged.update(r)
    tail = [ln for ln in (out + "\n" + err).splitlines() if not ln.startswith("@@RESULT@@")][-25:]
    return {"status": status, "rc": rc, "result": merged, "tail": tail}


def model_inis() -> list[Path]:
    return sorted(p for p in PARAMS.glob("*.ini") if "base_params" not in p.name)

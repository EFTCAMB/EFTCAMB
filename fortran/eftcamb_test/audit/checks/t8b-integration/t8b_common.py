"""Shared helpers for track t8b-integration checks: run a function in a subprocess importing camb from a tree."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
PARAMS = REPO / "fortran" / "eftcamb_test" / "parameters"
AUDIT = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
    "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
)
TREES = {n: AUDIT / "trees" / n for n in ("A0", "AD", "ADD", "B0", "B1")}
WORK = AUDIT / "work" / "t8b-integration"


def tree_path(t):
    return TREES.get(str(t), Path(t))


def run_in_tree(tree, script: Path, args: dict, threads: int = 1, timeout: float = 900) -> dict:
    """Run ``script`` (which must define main(args) -> dict) in a subprocess with camb imported from ``tree``."""
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    env.pop("PYTHONPATH", None)
    code = (
        "import sys, json; sys.path.insert(0, %r); sys.path.insert(1, %r);"
        "import runpy; mod = runpy.run_path(%r);"
        "res = mod['main'](json.loads(sys.argv[1])); print('@@RESULT@@' + json.dumps(res))"
    ) % (str(tree_path(tree)), str(HERE), str(script))
    WORK.joinpath("out").mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(HERE.parent))
    import harness  # machine-wide CPU slots (resource rules)

    with harness.slot(threads):
        p = subprocess.run(
            [sys.executable, "-c", code, json.dumps(args)],
            capture_output=True,
            text=True,
            env=env,
            timeout=timeout,
            cwd=WORK / "out",
            preexec_fn=lambda: os.nice(10),
        )
    for line in p.stdout.splitlines():
        if line.startswith("@@RESULT@@"):
            return json.loads(line[len("@@RESULT@@"):])
    return {"error": f"rc={p.returncode}", "stdout_tail": p.stdout[-3000:], "stderr_tail": p.stderr[-3000:]}


def read_model(camb, ini: str):
    """read an ini (absolute or relative to the test parameters dir) with the EFT keys (A0 tree read_ini)."""
    path = Path(ini)
    if not path.is_absolute():
        path = PARAMS / ini
    cp = camb.read_ini(str(path), no_validate=True)
    sys.path.insert(0, str(HERE.parent))
    import harness

    keys = harness.parse_ini(path)
    if hasattr(cp, "EFTCAMB") and cp.EFTCAMB.EFTflag == 0 and int(keys.get("EFTflag", "0")) != 0:
        # pre-rebase tree (B1): read_ini ignores EFTCAMB keys
        cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)
    return cp

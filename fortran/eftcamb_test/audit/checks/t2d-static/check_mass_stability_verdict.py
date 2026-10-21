"""Stability verdict with the mass-stability criterion, code under audit (A0) vs a tree with
compute_additional_derivs fixed (FIX1).

For every test model (and the audit inis of models without a test ini) the background + stability check is run
(camb.get_background -> CAMBdata_SetParams -> EFTCAMB_Stability_Check) with EFT_mass_stability = T forced, in a fresh
process per (tree, model).  Prints the verdict on each tree and lists the models whose verdict differs.

Usage: python check_mass_stability_verdict.py OUT.json [TREE_A TREE_B] [--only-own]
  --only-own: do not force EFT_mass_stability, use each ini's own setting.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
PAR = REPO / "fortran/eftcamb_test/parameters"
AUD = REPO / "fortran/eftcamb_test/audit/inis"
SCR = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit")
CHILD = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
from camb.baseconfig import CAMBError
pars = camb.read_ini(sys.argv[2], no_validate=True)
try:
    camb.get_background(pars)
    print("JSON" + json.dumps({"verdict": "stable"}))
except CAMBError as e:
    msg = str(e).strip().splitlines()[-1]
    print("JSON" + json.dumps({"verdict": "unstable" if "unstable" in msg else "error", "msg": msg}))
"""


def run(tree: Path, ini: Path, tmp: Path, force: bool) -> dict:
    w = tmp / f"{tree.name}_{ini.stem}.ini"
    lines = [f"DEFAULT({ini})", "feedback_level = 0", "EFTCAMB_feedback_level = 0", "EFTCAMB_write_background = F",
             f"output_root = {tmp / (tree.name + '_' + ini.stem)}"]
    if force:
        lines.append("EFT_mass_stability = T")
    w.write_text("\n".join(lines) + "\n")
    try:
        r = subprocess.run([sys.executable, "-c", CHILD, str(tree), str(w)], capture_output=True, text=True, cwd=tmp,
                           env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900)
    except subprocess.TimeoutExpired:
        return {"verdict": "timeout"}
    js = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
    if not js:
        return {"verdict": f"crash rc={r.returncode}", "tail": (r.stdout + r.stderr)[-300:]}
    return json.loads(js[-1][4:])


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    force = "--only-own" not in sys.argv
    out = Path(args[0])
    trees = [Path(t) for t in args[1:3]] or [SCR / "trees/A0", SCR / "work/t2d-static/FIX1"]
    inis = sorted(PAR.glob("[0-9]*.ini")) + sorted(AUD.glob("00[78]p*.ini"))
    tmp = Path(tempfile.mkdtemp(prefix="massstab_"))
    jobs = [(t, i) for i in inis for t in trees]
    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(lambda ti: run(ti[0], ti[1], tmp, force), jobs))
    table = {}
    for (t, i), r in zip(jobs, res):
        table.setdefault(i.stem, {})[t.name] = r
    ndiff = 0
    for m, d in table.items():
        v = [d[t.name]["verdict"] for t in trees]
        flag = "  <== DIFFERENT" if len(set(v)) > 1 else ""
        ndiff += bool(flag)
        print(f"{m:34s} " + "  ".join(f"{t.name}={x:10s}" for t, x in zip(trees, v)) + flag)
    print(f"models: {len(table)}, verdict differs: {ndiff}  (force EFT_mass_stability = {force})")
    out.write_text(json.dumps(table, indent=1))


if __name__ == "__main__":
    main()

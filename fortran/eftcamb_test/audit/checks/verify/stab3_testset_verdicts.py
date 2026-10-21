"""Does the loop-3 y fix change the stability verdict of any shipped test model?

For every fortran/eftcamb_test/parameters/*.ini (73 models) run camb.get_background (stability check on, as in the
default pipeline) on A0 and on FIXY (A0 + the one-line y update) and compare the verdicts / error messages.

Usage: python stab3_testset_verdicts.py [chunk_index n_chunks]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45"
             "/scratchpad/audit")
TREES = {"A0": AUDIT / "trees/A0", "FIXY": AUDIT / "work/t6-stability/fixy"}
OUT = AUDIT / "work/verify/stab3"

WORKER = r"""
import sys, json, time
sys.path.insert(0, sys.argv[1])
import camb
for w in sys.argv[2:]:
    t = time.time()
    try:
        cp = camb.read_ini(w, no_validate=True)
        camb.get_background(cp)
        v = "stable"
    except Exception as e:
        v = "error: " + str(e).strip().splitlines()[-1][:200]
    print("RESULT " + json.dumps({"ini": w, "verdict": v, "t": time.time() - t}), flush=True)
"""


def main():
    inis = harness.model_inis()
    if len(sys.argv) > 2:
        ci, nc = int(sys.argv[1]), int(sys.argv[2])
        inis = inis[ci::nc]
    OUT.mkdir(parents=True, exist_ok=True)
    store = OUT / "testset_verdicts.json"
    allres = json.loads(store.read_text()) if store.exists() else {}
    with tempfile.TemporaryDirectory(dir=OUT) as d:
        wrapped = []
        for ini in inis:
            w = Path(d) / Path(ini).name
            w.write_text(f"DEFAULT({ini})\nfeedback_level = 0\nEFTCAMB_write_background = F\n"
                         f"output_root = {Path(d) / Path(ini).stem}\n")
            wrapped.append(str(w))
        for tname, tree in TREES.items():
            env = dict(os.environ, OMP_NUM_THREADS="1")
            env.pop("PYTHONPATH", None)
            with harness.slot(1):
                p = subprocess.run([sys.executable, "-c", WORKER, str(tree), *wrapped], capture_output=True,
                                   text=True, cwd=d, env=env, timeout=540, preexec_fn=lambda: os.nice(10))
            got = [json.loads(line[7:]) for line in p.stdout.splitlines() if line.startswith("RESULT ")]
            if len(got) != len(wrapped):
                print(tname, "incomplete:", len(got), "of", len(wrapped), p.returncode, p.stderr[-1500:])
            for g in got:
                allres.setdefault(Path(g["ini"]).stem, {})[tname] = g["verdict"]
    store.write_text(json.dumps(allres, indent=1, sort_keys=True))
    diff = {k: v for k, v in allres.items() if v.get("A0") != v.get("FIXY")}
    nonstable = {k: v for k, v in allres.items() if v.get("A0") != "stable"}
    print(f"models with both verdicts: {sum(1 for v in allres.values() if 'A0' in v and 'FIXY' in v)}")
    print("verdict differs A0 vs FIXY:", json.dumps(diff, indent=1))
    print("A0 non-stable verdicts:", json.dumps(nonstable, indent=1))


if __name__ == "__main__":
    main()

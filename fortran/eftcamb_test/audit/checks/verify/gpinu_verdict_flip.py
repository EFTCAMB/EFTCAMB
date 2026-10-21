"""Verifier: does the tripled gpinu_tot in the stability cache flip an end-to-end stability verdict?

For each model ini and each stability-flag configuration, runs camb.get_background with the stability check enabled
(EFTCAMB_skip_stability = F) in two trees and records stable / 'EFTCAMB: theory unstable' / other error.

Usage: python gpinu_verdict_flip.py OUT.json TREE_BUGGY TREE_FIXED model.ini [model.ini ...]
"""

import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

CONFIGS = {
    "ini": {},
    "mass_math": {"EFT_mass_math_stability": True},
    "mass": {"EFT_mass_stability": True},
    "positivity": {"EFT_positivity_bounds": True},
    "minkowski": {"EFT_minkowski_limit": True},
}

WORKER = r"""
import json, sys
sys.path.insert(0, sys.argv[1])
import camb
cfgs = json.loads(sys.argv[3])
out = {}
for name, flags in cfgs.items():
    cp = camb.read_ini(sys.argv[2], no_validate=True)
    cp.EFTCAMB.EFTCAMB_skip_stability = False
    cp.EFTCAMB.EFTCAMB_feedback_level = 0
    for k, v in flags.items():
        setattr(cp.EFTCAMB, k, v)
    try:
        camb.get_background(cp)
        out[name] = "stable"
    except Exception as e:  # noqa
        msg = str(e).strip().splitlines()[-1]
        out[name] = "unstable" if "unstable" in msg else "error: " + msg[:80]
print("RESULT " + json.dumps(out))
"""


def run(tree, ini):
    with tempfile.TemporaryDirectory() as d:
        w = Path(d) / "m.ini"
        w.write_text(f"DEFAULT({ini})\noutput_root = {d}/out\nfeedback_level = 0\nEFTCAMB_write_background = F\n")
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        with harness.slot(1):
            p = subprocess.run([sys.executable, "-c", WORKER, tree, str(w), json.dumps(CONFIGS)], capture_output=True,
                               text=True, cwd=d, env=env, preexec_fn=lambda: os.nice(10), timeout=900)
        res = [json.loads(line[7:]) for line in p.stdout.splitlines() if line.startswith("RESULT ")]
        if not res:
            return {"crash": (p.stdout[-600:] + p.stderr[-600:])}
        return res[0]


def main():
    out_json, tb, tf = sys.argv[1:4]
    inis = sys.argv[4:]
    jobs = [(t, i) for i in inis for t in (tb, tf)]
    with ThreadPoolExecutor(3) as ex:
        res = list(ex.map(lambda ti: run(*ti), jobs))
    table = {}
    nflip = 0
    for k, ini in enumerate(inis):
        rb, rf = res[2 * k], res[2 * k + 1]
        name = Path(ini).stem
        table[name] = {"buggy": rb, "fixed": rf}
        flips = [c for c in CONFIGS if rb.get(c) != rf.get(c)]
        nflip += len(flips)
        print(f"{name:32s} " + " ".join(f"{c}:{rb.get(c, '?')[:8]}/{rf.get(c, '?')[:8]}" for c in CONFIGS)
              + (f"   FLIP {flips}" if flips else ""))
    print("total verdict flips buggy->fixed:", nflip)
    Path(out_json).write_text(json.dumps(table, indent=1))


if __name__ == "__main__":
    main()

"""t3-gr (a): bit-identity test with FMA contraction disabled (-O3 -ffp-contract=off) in both builds.

B0nofma = plain CAMB v2; R1nofma = code under audit with max_redshift=50, original polterddot and Transfer_Weyl.
Also R1nofma vs A0-equivalent? (run with python_set max_redshift only, no source reverts) is covered by A0nofma if built.
python a8_nofma.py TREE_A TREE_B [case,...]
"""
import sys
from concurrent.futures import ThreadPoolExecutor
from common import harness, fmt, run, WORK
from a5_variants import CASES

ta, tb = sys.argv[1], sys.argv[2]
cases = sys.argv[3].split(",") if len(sys.argv) > 3 else ["base"]
sub = f"a8_{ta}_{tb}"
tree = {t: (t if t in harness.TREES else WORK / t) for t in (ta, tb)}
from common import GR_INI, RUNS


def run_cached(t, c):
    out = RUNS / sub / f"{t}_{c}.npz"
    if out.exists():
        return out
    ov = {k: v for k, v in CASES[c].items() if k != "__ini__"}
    return run(tree[t], f"{t}_{c}", ini=CASES[c].get("__ini__", GR_INI), overrides=ov, sub=sub,
               python_set={"Reion.max_redshift": 50.0})


with ThreadPoolExecutor(int(__import__("os").environ.get("T3_JOBS", "4"))) as ex:
    futs = {(c, t): ex.submit(run_cached, t, c) for c in cases for t in (ta, tb)}
    for c in cases:
        r = harness.compare(futs[(c, ta)].result(), futs[(c, tb)].result())
        print(f"{c:12s}", fmt(r))
        for k, d in r["cl"].items():
            print("   ", k, {s: f"{x:.1e}" for s, x in d.items()})

"""t3-gr (a): compare a variant tree with B0 on the GR ini (threads=2)."""
import sys
from common import harness, run, fmt, WORK

name = sys.argv[1]
sub = "a3"
tree = WORK / name
b0 = run("B0", "B0_GR", sub=sub)
v = run(tree, f"{name}_GR", sub=sub)
r = harness.compare(v, b0)
print(f"{name} vs B0:", fmt(r))
for k, d in r["cl"].items():
    print("  ", k, {s: f"{x:.1e}" for s, x in d.items()})

"""t3-gr (a): test whether the reionization max_redshift change explains the EFTflag=0 A0 vs B0 difference."""
from common import harness, run, fmt, RUNS
import sys

sub = "a1"
b0 = run("B0", "B0_GR", sub=sub)
a0 = run("A0", "A0_GR", sub=sub)
a0r = run("A0", "A0_GR_maxz50", python_set={"Reion.max_redshift": 50.0}, sub=sub)
b0r = run("B0", "B0_GR_maxz1000", python_set={"Reion.max_redshift": 1000.0}, sub=sub)
print("A0 vs B0           :", fmt(harness.compare(a0, b0)))
print("A0(maxz=50) vs B0  :", fmt(harness.compare(a0r, b0)))
print("A0 vs B0(maxz=1000):", fmt(harness.compare(a0, b0r)))
for n, p in [("B0", b0), ("A0", a0), ("A0r", a0r), ("B0r", b0r)]:
    d = harness.load(p)
    print(n, "zrei=", d["meta"].get("derived", {}).get("zrei"))

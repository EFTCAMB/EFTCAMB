"""t3-gr (a): run-to-run determinism of B0 and A0 (threads=2 and threads=1) to know whether bit-identity is testable."""
from concurrent.futures import ThreadPoolExecutor
from common import harness, run, fmt

sub = "a2"
ps = {"Reion.max_redshift": 50.0}
jobs = {
    "B0_t2_r1": ("B0", None, 2), "B0_t2_r2": ("B0", None, 2),
    "B0_t1": ("B0", None, 1), "A0_t1_maxz50": ("A0", ps, 1), "A0_t2_maxz50": ("A0", ps, 2),
}
with ThreadPoolExecutor(4) as ex:
    futs = {n: ex.submit(run, t, n, python_set=p, threads=th, sub=sub) for n, (t, p, th) in jobs.items()}
    out = {n: f.result() for n, f in futs.items()}
print("B0 t2 run1 vs run2 :", fmt(harness.compare(out["B0_t2_r1"], out["B0_t2_r2"])))
print("B0 t1 vs t2        :", fmt(harness.compare(out["B0_t1"], out["B0_t2_r1"])))
print("A0m t1 vs B0 t1    :", fmt(harness.compare(out["A0_t1_maxz50"], out["B0_t1"])))
print("A0m t2 vs B0 t2    :", fmt(harness.compare(out["A0_t2_maxz50"], out["B0_t2_r1"])))
r = harness.compare(out["A0_t1_maxz50"], out["B0_t1"])
for k, d in r["cl"].items():
    print("  ", k, {s: f"{v:.1e}" for s, v in d.items()})

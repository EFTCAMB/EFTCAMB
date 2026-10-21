"""Verifier (fullmap-bg-eft-stale-background): run BH / Ext. Galileon on the ITER tree with EFT_ITER = 0, 1, 5, 20
and compare spectra; parse VITMAX lines (used vs self-consistently iterated EFTc, EFTLambda, adotoa)."""

import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = Path(os.environ["AUDIT_ROOT"]) / "work" / "verify" / "stalebg"
TREE = W / "ITER"
model = sys.argv[1]  # bh | eg
ini = W / ("008p5-beyond-horndeski.ini" if model == "bh" else "008p7-ext-galileon.ini")
iters = [int(x) for x in sys.argv[2].split(",")]
trees = {f"iter{n}": (TREE, n) for n in iters}
if len(sys.argv) > 3:
    trees["A0"] = ("A0", 0)


def job(item):
    name, (tree, n) = item
    out = W / f"{model}_{name}.npz"
    st = harness.run_model(tree, ini, out, overrides={"EFTCAMB_skip_stability": "T"}, threads=1,
                           extra_env={"EFT_ITER": str(n)})
    return name, st, out


with ThreadPoolExecutor(4) as ex:
    res = list(ex.map(job, trees.items()))
for name, st, out in res:
    print(name, st["status"])
    vit = {}
    try:
        for line in Path(st["log"]).read_text().splitlines():
            if line.startswith("VITMAX"):
                v = line.split()
                vit[int(v[1])] = [float(x) for x in v[2:]]
    except Exception as e:  # noqa: BLE001
        print("  log read", e)
    for i, lab in ((1, "EFTc"), (2, "EFTLambda"), (3, "adotoa")):
        if i in vit:
            a, k, old, new, rel = vit[i]
            print(f"  {lab:10s} max rel(used vs iterated) {rel:.3e} at a={a:.4f} k={k:.3e}: {old:+.6e} -> {new:+.6e}")
    if st["status"] != "ok":
        print("  ", st.get("log_tail", "")[-800:])
ok = {n: o for n, s, o in res if s["status"] == "ok"}
names = list(ok)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        c = harness.compare(harness.load(ok[names[i]]), harness.load(ok[names[j]]))
        print(f"{names[i]} vs {names[j]}: worst_cl {c['worst_cl']} worst_pk {c['worst_pk']} "
              f"pp {json.dumps(c['cl'].get('lens_potential', ''), default=str)[:200]}")

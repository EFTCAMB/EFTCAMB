"""Run shipped models on the instrumented release tree (interp_xfinal_instr) and count out-of-bounds index hits.

mode 'scan'  : run models [start:end] as shipped (no clamp); report AUDIT_OOB hits (x, x(N), garbage x(N+1), y(N+1)).
mode 'clamp' : rerun the given models with EFT_AUDIT_CLAMP=1 and compare to the unclamped run of the same binary and to
               the frozen A0 reference runs.
Usage: python interp_xfinal_release_scan.py scan START END | clamp MODEL [MODEL ...]
"""

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness  # noqa: E402

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
TREE = AUDIT / "work/verify/interp_xfinal_instr"
OUT = AUDIT / "work/verify/interp_xfinal/rel_scan"
OUT.mkdir(parents=True, exist_ok=True)
inis = harness.model_inis()
by_name = {p.stem: p for p in inis}


def hits(log: Path) -> list[str]:
    return [ln for ln in log.read_text(errors="replace").splitlines() if ln.startswith("AUDIT_OOB")]


mode = sys.argv[1]
if mode == "scan":
    a, b = int(sys.argv[2]), int(sys.argv[3])
    todo = [p for p in inis[a:b] if not (OUT / f"{p.stem}.json").exists()]
    if todo:
        harness.run_batch(TREE, todo, OUT, jobs=4, threads=1, timeout=540)
    for p in inis:
        js = OUT / f"{p.stem}.json"
        if not js.exists():
            continue
        info = json.loads(js.read_text())
        h = hits(Path(info["log"]))
        line = f"{p.stem:36s} {info['status']:12s} oob_hits={len(h)}"
        if h:
            xs = [float(re.search(r" x=\s*(\S+)", s).group(1)) for s in h]
            g = sorted({re.search(r"x\(N\+1\)=\s*(\S+) y\(N\+1\)=\s*(\S+)", s).groups() for s in h})
            line += f" x in [{min(xs):.3e},{max(xs):.3e}] garbage(x(N+1),y(N+1)) distinct={len(g)} e.g. {g[:3]}"
        if info["status"] == "ok" and (AUDIT / "runs/A0" / f"{p.stem}.npz").exists():
            c = harness.compare(harness.load(OUT / f"{p.stem}.npz"), harness.load(AUDIT / "runs/A0" / f"{p.stem}.npz"))
            line += f" vsA0: bit_identical={c['bit_identical']} worst_cl={c['worst_cl']:.2e}"
        print(line)
elif mode == "clamp":
    names = sys.argv[2:]
    harness.run_batch(TREE, [by_name[n] for n in names], OUT, jobs=4, threads=1, timeout=540, suffix="_clamp",
                      extra_env={"EFT_AUDIT_CLAMP": "1"})
    for n in names:
        info = json.loads((OUT / f"{n}_clamp.json").read_text())
        h = hits(Path(info["log"]))
        line = f"{n:36s} clamp status={info['status']} oob_hits={len(h)}"
        if info["status"] == "ok" and (OUT / f"{n}.npz").exists():
            ca = harness.load(OUT / f"{n}_clamp.npz")
            cb = harness.load(OUT / f"{n}.npz")
            c = harness.compare(ca, cb)
            eft = max(c["eft"].values()) if c["eft"] else 0.0
            line += (f" clamp-vs-noclamp: bit_identical={c['bit_identical']} worst_cl={c['worst_cl']:.2e} "
                     f"worst_pk={c['worst_pk']:.2e} sigma8={c.get('sigma8', float('nan')):.2e} max_eft={eft:.2e}")
        print(line)

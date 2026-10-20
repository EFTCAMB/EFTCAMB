"""Non-GR follow-up for t3-3p1-success-unset: do the 3+1 test inis run once the background flag is fixed?

The test inis set EFT_ghost_stability=T, which EFTCAMB rejects for beyond-GLPV models by design ('stability for model
beyond GLPV has not been worked out ... disable EFT_ghost_stability'). Here the generic ghost/gradient checks are off,
the model's own Designer3p1Ghost/GradientStability checks stay on (defaults T).
Runs A0, B1, R7 on both inis; compares R7 vs B1 (CAMB 1.6.5 -> v2 shape expected) and R7 vs A0 EFTflag=0 (size of EFT
effect) at full default settings.
Usage: python verify_3p1_nonGR.py R7_TREE OUTDIR
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

R7 = sys.argv[1]
OUT = Path(sys.argv[2])
INIS = REPO / "fortran" / "eftcamb_test" / "audit" / "inis"
OV = {"EFT_ghost_stability": "F", "EFT_gradient_stability": "F"}
JOBS = [(t, f"{tag}_{n}_nostab", INIS / f"007p8-3p1{s}.ini") for n, s in (("ini1", ""), ("ini2", "_2"))
        for t, tag in (("A0", "A0"), ("B1", "B1"), (R7, "R7"))]


def job(args):
    tree, name, ini = args
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, ini, out, overrides=OV, threads=1, timeout=900)
    log = out.with_suffix(".log").read_text(errors="replace").splitlines()
    msg = [ln.strip() for ln in log if "EFTCAMB:" in ln or "DLSODA" in ln or "Instability" in ln][-2:]
    return name, st, msg


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(4) as ex:
        res = {n: (st, msg) for n, st, msg in ex.map(job, JOBS)}
    summary = {}
    for n, (st, msg) in res.items():
        summary[n] = {"status": st["status"], "msg": msg}
        print(f"{n:22s} {st['status']:12s} {msg}")
    lcdm = OUT.parent / "runs" / "A0_LCDM.npz"
    for n in ("ini1", "ini2"):
        r7, b1 = OUT / f"R7_{n}_nostab.npz", OUT / f"B1_{n}_nostab.npz"
        if r7.exists() and b1.exists():
            c = harness.compare(r7, b1)
            summary[f"R7_vs_B1_{n}"] = c["cl"]
            print(f"R7 vs B1 {n}: worst_cl={c['worst_cl']:.2e} worst_pk={c['worst_pk']:.2e} "
                  f"TT={c['cl']['total']['TT']:.2e} EE={c['cl']['total']['EE']:.2e} "
                  f"tensorBB={c['cl'].get('tensor', {}).get('BB', float('nan')):.2e} sigma8={c.get('sigma8', -1):.2e}")
        if r7.exists() and lcdm.exists():
            c = harness.compare(r7, lcdm)
            print(f"R7 {n} vs A0 EFTflag=0: TT={c['cl']['total']['TT']:.2e} sigma8={c.get('sigma8', -1):.2e}")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()

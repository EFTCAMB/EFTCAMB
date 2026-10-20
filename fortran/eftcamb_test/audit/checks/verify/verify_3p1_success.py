"""Independent check of finding t3-3p1-success-unset (3+1 wDE designer always 'background solver failed').

Runs, at the DEFAULT (full) settings of the inis (lensing, tensors, l_max as in base_params), not the FAST settings of
the repro:
  - A0 and B1 on both 3+1 test inis (007p8-3p1.ini, 007p8-3p1_2.ini)
  - the fixed tree R7 (only change: success = .true. at top of SolveDesignerEquations) on both inis
  - R7 at the GR point (sigma1 = rho1 = 1, w = -1) with the EFT path forced and stability skipped, vs A0 EFTflag=0
    (1_EFT_GR.ini) at IntTolBoost = 1 and 2.
Usage: python verify_3p1_success.py R7_TREE OUTDIR
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
I1, I2 = INIS / "007p8-3p1.ini", INIS / "007p8-3p1_2.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
GRPT = {"EFTwDE": "0", "Designer3p1modelsigma1": "1", "Designer3p1sigma10": "1.0", "Designer3p1modelsigma1_ODE": "0",
        "Designer3p1modelrho1N": "0", "Designer3p1modelrho2": "0", "Designer3p1LuminalGW": "T"}
FORCED = {"EFTCAMB_GR_threshold": "-1", "EFTCAMB_skip_stability": "T"}

JOBS = [
    ("A0", "A0_ini1", I1, {}, {}),
    ("A0", "A0_ini2", I2, {}, {}),
    ("B1", "B1_ini1", I1, {}, {}),
    ("B1", "B1_ini2", I2, {}, {}),
    (R7, "R7_ini1", I1, {}, {}),
    (R7, "R7_ini2", I2, {}, {}),
    (R7, "R7_GRforced", I1, {**GRPT, **FORCED}, {}),
    (R7, "R7_GRforced_tol2", I1, {**GRPT, **FORCED}, {"Accuracy.IntTolBoost": 2.0}),
    ("A0", "A0_LCDM", GR, {}, {}),
    ("A0", "A0_LCDM_tol2", GR, {}, {"Accuracy.IntTolBoost": 2.0}),
]


def job(args):
    tree, name, ini, ov, ps = args
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, ini, out, overrides=ov, python_set=ps, threads=1, timeout=900)
    log = out.with_suffix(".log").read_text(errors="replace").splitlines()
    msg = [ln.strip() for ln in log if "EFTCAMB:" in ln or "DLSODA" in ln or "unstable" in ln.lower()][-2:]
    return name, st, msg


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(4) as ex:
        res = {n: (st, msg) for n, st, msg in ex.map(job, JOBS)}
    summary = {}
    for n, (st, msg) in res.items():
        row = {"status": st["status"], "msg": msg}
        if st["status"] == "ok" and n.startswith("R7_GRforced"):
            ref = OUT / ("A0_LCDM_tol2.npz" if n.endswith("tol2") else "A0_LCDM.npz")
            c = harness.compare(OUT / f"{n}.npz", ref)
            row.update(worst_cl=c["worst_cl"], worst_pk=c["worst_pk"], sigma8=c["sigma8"], ref=ref.name,
                       cl={k: v for k, v in c["cl"].items()})
        if st["status"] == "ok":
            meta = harness.load(OUT / f"{n}.npz").get("meta", {})
            row["eft_model_name"] = meta.get("eft_model_name") if isinstance(meta, dict) else None
        summary[n] = row
        extra = f"  worst_cl={row['worst_cl']:.2e} worst_pk={row['worst_pk']:.2e}" if "worst_cl" in row else ""
        print(f"{n:20s} {st['status']:12s} {msg}{extra}")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()

"""K-mimic (5_Kmimic_1, 5_Kmimic_3): spectra vs model_background_num_points and vs small accuracy changes.
Tests whether the large default-vs-10x change seen in bgtable_testmodels.py is a table-resolution effect or a
hypersensitivity of the K-mimic perturbations to 1e-4 background changes (which would make the proposed fix
'raise num_points to 10000' unsafe for K-mimic)."""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

P = harness.PARAMS_DIR
OUT = harness.AUDIT_ROOT / "work" / "verify" / "bgtable_kmimic"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "2000", "num_redshiftwindows": "0"}
JOBS = {}
for m in ("5_Kmimic_1", "5_Kmimic_3"):
    for n in (1000, 1001, 1500, 2000, 3000, 4000, 6000, 10000, 20000):
        JOBS[f"{m}__N{n}"] = (m, {"model_background_num_points": str(n)}, {})
    JOBS[f"{m}__N1000_itb2"] = (m, {}, {"Accuracy.IntTolBoost": 2.0})
    JOBS[f"{m}__N10000_itb2"] = (m, {"model_background_num_points": "10000"}, {"Accuracy.IntTolBoost": 2.0})


def job(item):
    name, (m, ov, ps) = item
    out = OUT / f"{name}.npz"
    if out.exists() and json.loads(out.with_suffix(".json").read_text())["status"] == "ok":
        return name, {"status": "ok"}
    return name, harness.run_model("A0", P / f"{m}.ini", out, overrides={**FAST, **ov}, python_set=ps, threads=1,
                                   timeout=1500)


def main() -> None:
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(job, JOBS.items()))
    print(f"{'run':28s} {'status':8s} {'sigma8':>9s} {'TT200':>10s} {'TT1000':>10s} {'EE500':>10s}")
    for n in JOBS:
        st = res[n]
        if st["status"] != "ok":
            print(f"{n:28s} {st['status']}")
            continue
        d = harness.load(OUT / f"{n}.npz")
        c = d["cl_unlensed_scalar"]
        print(f"{n:28s} {'ok':8s} {float(d['sigma8'][0]):9.5f} {c[200, 0]:10.4e} {c[1000, 0]:10.4e} {c[500, 1]:10.4e}")
    for m in ("5_Kmimic_1", "5_Kmimic_3"):
        d = harness.load(harness.AUDIT_ROOT / "runs" / "A0" / f"{m}.npz")
        c = d["cl_unlensed_scalar"]
        print(f"{m + ' A0 ref (full ini)':28s} {'':8s} {float(d['sigma8'][0]):9.5f} {c[200, 0]:10.4e} "
              f"{c[1000, 0]:10.4e} {c[500, 1]:10.4e}")
        b1 = harness.AUDIT_ROOT / "runs" / "B1" / f"{m}.npz"
        if b1.exists():
            d = harness.load(b1)
            c = d["cl_unlensed_scalar"]
            print(f"{m + ' B1 ref (full ini)':28s} {'':8s} {float(d['sigma8'][0]):9.5f} {c[200, 0]:10.4e} "
                  f"{c[1000, 0]:10.4e} {c[500, 1]:10.4e}")
    _ = np


if __name__ == "__main__":
    main()

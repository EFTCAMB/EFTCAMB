"""Self-convergence of the background-table resolution for the shipped test models (not the GR limit):
default model_background_num_points vs 10x, same model and settings. Tells whether the table error in genuinely
modified models is larger or smaller than the ~1e-4 GR-limit floor reported in t3-bg-table-resolution-gr-floor.
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

P = harness.PARAMS_DIR
OUT = harness.AUDIT_ROOT / "work" / "verify" / "bgtable_models"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "2000", "num_redshiftwindows": "0"}
ADE_SET = {"get_tensor_cls": "F"}
MODELS = {
    "5_Kmouflage_1": 1000, "5_Kmouflage_2": 1000, "5_Kmouflage_3": 1000,
    "5_Kmimic_1": 1000, "5_Kmimic_3": 1000,
    "5_ADE_1": 1500,
    "5_quint_axion": 1100, "5_quint_cubic": 1100, "5_quint_exponential": 1100, "5_quint_inverse_power": 1100,
    "5_quint_quadratic": 1100,
}


def job(item):
    name, ov = item
    model = name.split("__")[0]
    base = ADE_SET if "ADE" in model else FAST
    out = OUT / f"{name}.npz"
    if out.exists() and json.loads(out.with_suffix(".json").read_text())["status"] == "ok":
        return name, {"status": "ok"}
    return name, harness.run_model("A0", P / f"{model}.ini", out, overrides={**base, **ov}, threads=1, timeout=1500)


def main() -> None:
    jobs = []
    for m, n in MODELS.items():
        jobs.append((f"{m}__N{n}", {}))
        jobs.append((f"{m}__N{10 * n}", {"model_background_num_points": str(10 * n)}))
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(job, jobs))
    for n, st in res.items():
        if st["status"] != "ok":
            print(n, st["status"], st.get("log_tail", [])[-5:])
    rows = []
    print(f"{'model':24s} {'dH/H':>9s} {'TT':>9s} {'EE':>9s} {'TE':>9s} {'P(k)':>9s}   (default N vs 10x N)")
    for m, n in MODELS.items():
        a, b = OUT / f"{m}__N{n}.npz", OUT / f"{m}__N{10 * n}.npz"
        if not (a.exists() and b.exists()):
            continue
        r = harness.compare(a, b, lmax=1900)
        u = r["cl"]["unlensed_scalar"]
        rows.append(dict(model=m, bg_H=r["bg_H"], TT=u["TT"], EE=u["EE"], TE=u["TE"], pk=r["worst_pk"]))
        print(f"{m:24s} {r['bg_H']:9.2e} {u['TT']:9.2e} {u['EE']:9.2e} {u['TE']:9.2e} {r['worst_pk']:9.2e}")
    (OUT / "summary.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()

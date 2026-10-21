"""Verifier for t3-3p1-success-unset: run 3+1 designer on A0, AD-like trees, B1 and fixed R7 copy."""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

INI = HERE / "v_3p1.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "v_3p1"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "200", "get_transfer": "F",
        "num_redshiftwindows": "0"}
GRPT = {"EFTwDE": "0", "Designer3p1modelsigma1": "1", "Designer3p1sigma10": "1.0", "Designer3p1modelsigma1_ODE": "0",
        "Designer3p1modelrho1N": "0", "Designer3p1modelrho2": "0"}
FORCED = {"EFTCAMB_GR_threshold": "-1", "EFTCAMB_skip_stability": "T"}
trees = sys.argv[1:] or ["A0"]


def job(a):
    tree, name, ini, ov, ps = a
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, ini, out, overrides={**FAST, **ov}, python_set=ps, threads=1, timeout=600)
    txt = out.with_suffix(".log").read_text(errors="replace")
    msg = [ln.strip() for ln in txt.splitlines() if "solver failed" in ln or "unstable" in ln or "GLPV" in ln
           or "DLSODA" in ln or "rror" in ln][-2:]
    return name, st, out, msg


jobs = [("A0", "A0_LCDM", GR, {}, {"Accuracy.IntTolBoost": 2.0})]
for t in trees:
    tn = Path(t).name
    jobs += [(t, f"{tn}_ini", INI, {}, {}), (t, f"{tn}_GR", INI, GRPT, {}),
             (t, f"{tn}_GRforced", INI, {**GRPT, **FORCED}, {})]
with ThreadPoolExecutor(3) as ex:
    res = {n: (st, out, msg) for n, st, out, msg in ex.map(job, jobs)}
ref = res["A0_LCDM"][1]
for n, (st, out, msg) in res.items():
    line = f"{n:22s} {st['status']:12s} {msg}"
    if st["status"] == "ok" and n != "A0_LCDM":
        line += f"  worst_cl vs LCDM: {harness.compare(out, ref)['worst_cl']:.2e}"
    print(line)

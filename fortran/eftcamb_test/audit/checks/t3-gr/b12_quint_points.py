"""t3-gr (b): quintessence (exponential potential) GR approach saturates at ~5e-5 in H for EFT_lambda <= 0.01 (b11).
Same background-table-resolution floor as K-mouflage/ADE? Rerun with 10x model_background_num_points (1100 -> 11000).
"""
from concurrent.futures import ThreadPoolExecutor

from common import PARAMS, RUNS, fmt, harness

SUB = RUNS / "b11"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "1500", "num_redshiftwindows": "0"}
INI = PARAMS / "5_quint_exponential.ini"
JOBS = {f"quint_exp_lam{lam:g}{'_np11000' if np else ''}": {"EFT_lambda": repr(lam),
                                                            **({"model_background_num_points": "11000"} if np else {})}
        for lam in (0.03, 0.01, 0.003) for np in (False, True)}


def one(item):
    n, ov = item
    out = SUB / f"{n}.npz"
    if not out.with_suffix(".json").exists():
        harness.run_model("A0", INI, out, overrides={**FAST, **ov}, threads=1, timeout=540)
    return n, out


with ThreadPoolExecutor(4) as ex:
    res = dict(ex.map(one, JOBS.items()))
ref = RUNS / "b10" / "REF.npz"
for n, out in res.items():
    st = harness.json.loads(out.with_suffix(".json").read_text())
    print(f"{n:28s} {st['status']:8s}", fmt(harness.compare(out, ref, lmax=1450)) if st["status"] == "ok" else "")

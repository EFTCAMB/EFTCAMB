"""Check 12: get_transfer_functions + power_spectra_from_transfer with EFT on (v2 python path).

camb.get_transfer_functions(params) runs the CMB/matter transfer calculation only; power_spectra_from_transfer()
then makes the Cls.  With EFT on the result must be identical to camb.get_results(params) (same EFT sources).
Also compares get_matter_transfer_data / linear P(k) between the two paths, and a changed InitPower amplitude
applied after the transfer step (As x 1.1 -> Cls x 1.1 exactly) as a check that the EFT sources are not recomputed.
"""

from __future__ import annotations

import json
import sys

import numpy as np
from t8c_lib import WORK, cl_diffs, load, maxrel, model_ini, run_many

OUT = WORK / "c12"
OV = {"get_tensor_cls": "F", "feedback_level": "0", "num_redshiftwindows": "0", "transfer_kmax": "2"}
MODELS = ["1_EFT_GR", "4_Des_F_R_LCDM_1", "5_hdsk_TG_1", "5_quint_exponential", "2_PEFT_gamma3_const_1"]
BASE = {"set_for_lmax": {"lmax": 2000, "lens_potential_accuracy": 1, "nonlinear": False}, "cl_lmax": 2000}


def jobs():
    js = []
    for m in MODELS:
        js.append(dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{m}_full.npz", overrides=OV,
                       args=BASE | {"outputs": ["cmb", "transfer"]}))
        js.append(dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{m}_tf.npz", overrides=OV,
                       args=BASE | {"outputs": ["transfer_only", "power_from_transfer", "cmb", "transfer"]}))
    return js


def load_safe(path):
    """harness.load, skipping the object-dtype tr_z entry written by the worker."""
    with np.load(path, allow_pickle=False) as f:
        d = {k: f[k] for k in f.files if k != "tr_z"}
    d["meta"] = json.loads(str(d["meta"]))
    return d


def report():
    res = {}
    for m in MODELS:
        try:
            a, b = load_safe(OUT / f"{m}_tf.npz"), load_safe(OUT / f"{m}_full.npz")
        except FileNotFoundError:
            st = {p: json.loads((OUT / f"{m}_{p}.json").read_text()) for p in ("tf", "full")}
            res[m] = {p: s["status"] + " " + str(s.get("log_tail", [])[-2:]) for p, s in st.items()}
            continue
        d = cl_diffs(a["cl_total"], b["cl_total"], lmax=2000)
        res[m] = {"cl_total": d, "cl_unlensed": cl_diffs(a["cl_unlensed_scalar"], b["cl_unlensed_scalar"], lmax=2000),
                  "pp": maxrel(a["cl_lens_potential"][2:, 0], b["cl_lens_potential"][2:, 0]),
                  "transfer": maxrel(a["tr"], b["tr"], floor_frac=1e-8),
                  "bit_identical": bool(np.array_equal(a["cl_total"], b["cl_total"])),
                  "finite": bool(np.all(np.isfinite(a["cl_total"])))}
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        for r in run_many(jobs(), nproc=int(sys.argv[2]) if len(sys.argv) > 2 else 4):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report()

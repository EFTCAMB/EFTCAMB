"""Release-build impact of the uninitialised grhog_t/grhor_t read in the EFT block of derivst (tensor equations).

Tree P = A0 + inert audit patches (make_patched_tree.py). AUDIT_FIX_TENSOR=1 assigns grhog_t, grhor_t before the
EFT tensor block. Compares tensor Cls (and scalar Cls as a control: must be bit-identical) for full-mapping models,
whose compute_H_derivs uses gpresdotdotm_t.

Usage: python tensor_fix.py
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/t8a-equations")
OUT = W / "tensor_fix"
OUT.mkdir(exist_ok=True)
MODELS = ["5_Horava_eta_4", "5_ADE_1", "5_quint_exponential"]
OV = {"l_max_scalar": "1500", "do_lensing": "F", "get_transfer": "F", "num_redshiftwindows": "0"}

jobs = [(m, tag, env) for m in MODELS for tag, env in (("off", {}), ("fix", {"AUDIT_FIX_TENSOR": "1"}))]


def run(j):
    m, tag, env = j
    return harness.run_model(W / "P", harness.PARAMS_DIR / f"{m}.ini", OUT / f"{m}_{tag}.npz", overrides=OV,
                             extra_env=env, threads=1)


with ThreadPoolExecutor(4) as ex:
    st = list(ex.map(run, jobs))
res = {}
for m in MODELS:
    a, b = OUT / f"{m}_off.npz", OUT / f"{m}_fix.npz"
    if not (a.exists() and b.exists()):
        res[m] = [s["status"] for s in st]
        continue
    A, B = harness.load(a), harness.load(b)
    ta, tb = A["cl_tensor"][2:1000], B["cl_tensor"][2:1000]
    with np.errstate(all="ignore"):
        r = {s: float(np.nanmax(np.abs(tb[:, i] - ta[:, i]) / np.abs(tb[:, i]))) for i, s in enumerate("TT EE BB".split())}
    r["scalar_bit_identical"] = bool(np.array_equal(A["cl_unlensed_scalar"], B["cl_unlensed_scalar"]))
    res[m] = r
    print(m, json.dumps(r))
(OUT / "tensor_fix.json").write_text(json.dumps(res, indent=1))

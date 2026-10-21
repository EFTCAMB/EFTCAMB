"""t3-gr (b): f(R) just below the RGR switch: B0 = 3e-8 is treated as pure GR by default (RGR time > 1). What does the
model give with the RGR switch disabled (EFTCAMB_GR_threshold = -1)? Matter power only (transfer_kmax = 10/Mpc, no
CMB) to keep the stiff (Compton mass ~ H/sqrt(B0)) pi evolution affordable.
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np

from common import GR_INI, PARAMS, RUNS, harness

INI = PARAMS / "4_Des_F_R_LCDM_1.ini"
SUB = RUNS / "b6"
PK = {"get_scalar_cls": "F", "get_tensor_cls": "F", "do_lensing": "F", "transfer_kmax": "10",
      "num_redshiftwindows": "0"}
JOBS = {
    "pk_GR": (GR_INI, PK, {"Accuracy.IntTolBoost": 2.0}),
    "pk_fR_default_1em07": (INI, {**PK, "EFTB0": "1e-07"}, {}),
    "pk_fR_default_3em08": (INI, {**PK, "EFTB0": "3e-08"}, {}),
    "pk_fR_forced_3em08": (INI, {**PK, "EFTB0": "3e-08", "EFTCAMB_GR_threshold": "-1"}, {}),
}


def one(item):
    n, (ini, ov, ps) = item
    out = SUB / f"{n}.npz"
    if not out.with_suffix(".json").exists():
        harness.run_model("A0", ini, out, overrides=ov, python_set=ps, threads=1, timeout=480)
    return n, out


with ThreadPoolExecutor(4) as ex:
    res = dict(ex.map(one, JOBS.items()))
ref = harness.load(res["pk_GR"])
k = ref["pk_k"]
ik = [int(np.argmin(abs(k - x))) for x in (0.1, 0.3, 1.0, 3.0, 7.0)]
for n, out in res.items():
    st = harness.json.loads(out.with_suffix(".json").read_text())
    if st["status"] != "ok":
        print(f"{n:22s} {st['status']} elapsed {st.get('elapsed')}")
        continue
    d = harness.load(out)
    r = d["pk"][0] / ref["pk"][0] - 1
    print(f"{n:22s} " + " ".join(f"P(k={k[i]:.2g})-1={r[i]:+.2e}" for i in ik)
          + f" sigma8-1={d['sigma8'][0] / ref['sigma8'][0] - 1:+.1e} elapsed {st['elapsed']}s")

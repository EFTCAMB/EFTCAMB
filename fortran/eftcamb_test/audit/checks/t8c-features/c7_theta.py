"""Check 7: H0 from thetastar / cosmomc_theta with EFT on (new in CAMB v2: ini keys thetastar, cosmomc_theta).

For each model and theta kind:
  ini path   : flattened model ini with 'hubble' removed and thetastar (or cosmomc_theta) set -> camb.read_ini
  python path: read_ini(model ini) then CAMBparams.set_H0_for_theta(theta)  (EFTCAMB already initialised)
Then the theta of the final results is compared with the requested value (brentq rtol on H0 is 5e-5).

usage: python c7_theta.py [run|report]
"""

from __future__ import annotations

import json
import sys

from t8c_lib import WORK, harness, load, model_ini, run_many

OUT = WORK / "c7"
MODELS = ["1_EFT_GR", "2_PEFT_wDE_wCDM_1", "4_Des_F_R_wCDM_1", "4_Des_F_R_LCDM_1", "5_quint_exponential",
          "5_Kmouflage_1", "5_hdsk_TG_1"]
TARGET = {"thetastar": 0.0104110, "cosmomc_theta": 0.0104090}
OV = {"get_tensor_cls": "F", "get_transfer": "F", "get_scalar_cls": "F", "feedback_level": "0", "num_redshiftwindows": "0"}


def theta_ini(model: str, kind: str) -> str:
    keys = harness.parse_ini(model_ini(model))
    keys.pop("hubble", None)
    keys |= OV
    keys[kind] = repr(TARGET[kind])
    keys["EFTCAMB_write_background"] = "F"
    path = OUT / "inis" / f"{model}_{kind}.ini"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(f"{k} = {v}" for k, v in keys.items()) + "\n")
    return str(path)


def jobs():
    js = []
    for m in MODELS:
        for kind in TARGET:
            js.append(dict(tree="A0", ini=theta_ini(m, kind), task="generic", out=OUT / f"{m}_{kind}_ini.npz",
                           args={"outputs": ["background_only", "theta"]}))
            js.append(dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{m}_{kind}_py.npz", overrides=OV,
                           args={"theta_via_set_cosmology": [kind, TARGET[kind]],
                                 "outputs": ["background_only", "theta"]}))
            js.append(dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{m}_{kind}_pydirect.npz",
                           overrides=OV, args={"set_H0_for_theta": {"theta": TARGET[kind],
                                                                    "cosmomc_approx": kind == "cosmomc_theta"},
                                               "outputs": ["background_only", "theta"]}))
    return js


def report():
    res = {}
    for m in MODELS:
        for kind in TARGET:
            row = {}
            for path in ("ini", "py", "pydirect"):
                try:
                    d = load(OUT / f"{m}_{kind}_{path}.npz")["meta"]
                except FileNotFoundError:
                    st = json.loads((OUT / f"{m}_{kind}_{path}.json").read_text())
                    row[path] = st["status"] + ": " + (st.get("log_tail") or [""])[-1][-80:]
                    continue
                got = d["thetastar"] if kind == "thetastar" else d["cosmomc_theta"]
                row[path] = {"H0": d["H0_final"], "theta": got, "rel_err": got / TARGET[kind] - 1}
            res[f"{m}:{kind}"] = row
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        for r in run_many(jobs(), nproc=int(sys.argv[2]) if len(sys.argv) > 2 else 4):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report()

"""Independent verification of finding ini-theta-h0-solved-without-eft (read_ini / CAMB_ReadParamFile path).

For each model:
  H0_ini   : H0 produced by camb.read_ini with thetastar (or cosmomc_theta) in the ini
  th_ini   : theta of the final results for those params
  H0_gr    : H0 produced by the same ini with EFTflag forced to 0 (LCDM / CAMB-DE background)
  th_hub   : theta from an ini with hubble = H0_ini (EFT active) -> consistency of th_ini
  H0_true  : own scipy brentq on hubble (each trial re-reads the ini, full EFT init), rtol 1e-9
Control: CAMB native fluid DE (EFTflag=0, w=-0.8) solved through ini thetastar.
Usage: python initheta_indep.py [model ...]
"""

import json
import os
import sys
import tempfile
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
REPO = Path("/Users/marco/IW/eftcamb_update/EFTCAMB")
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

TREE = Path(os.environ.get("INITHETA_TREE", "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
                                                "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0"))
OUT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit/work/verify/initheta")
OUT.mkdir(parents=True, exist_ok=True)
PARAMS = REPO / "fortran/eftcamb_test/parameters"
sys.path.insert(0, str(TREE))
import camb  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

assert Path(camb.__file__).resolve().is_relative_to(TREE), camb.__file__

THETA = 0.0104110
TMP = Path(tempfile.mkdtemp(dir=OUT))
QUIET = {"get_scalar_cls": "F", "get_tensor_cls": "F", "get_transfer": "F", "num_redshiftwindows": "0",
         "feedback_level": "0", "EFTCAMB_feedback_level": "0", "EFTCAMB_write_background": "F"}


def flatten(path: Path) -> dict:
    keys, defaults = {}, []
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.upper().startswith("DEFAULT("):
            defaults.append(flatten(path.parent / line[8:line.rindex(")")].strip()))
        elif "=" in line:
            k, v = line.split("=", 1)
            keys[k.strip()] = v.strip()
    merged = {}
    for d in defaults:
        merged.update(d)
    return merged | keys


_n = 0


def params_from(keys: dict):
    global _n
    _n += 1
    ini = TMP / f"p{_n}.ini"
    ini.write_text("\n".join(f"{k} = {v}" for k, v in keys.items()) + "\n")
    return camb.read_ini(str(ini))


def theta_of(cp, kind):
    res = camb.get_background(cp)
    if kind == "cosmomc_theta":
        return res.cosmomc_theta()
    return res.get_derived_params()["thetastar"] / 100


def check(model: str, kind: str = "thetastar", extra: dict | None = None) -> dict:
    base = flatten(PARAMS / f"{model}.ini") | QUIET | (extra or {})
    base.pop("hubble", None)
    th_keys = base | {kind: str(THETA)}
    cp = params_from(th_keys)
    H0_ini = cp.H0
    th_ini = theta_of(cp, kind)
    cp_gr = params_from(th_keys | {"EFTflag": "0"})
    H0_gr = cp_gr.H0
    th_hub = theta_of(params_from(base | {"hubble": repr(H0_ini)}), kind)

    def f(h):
        return theta_of(params_from(base | {"hubble": repr(h)}), kind) - THETA

    # bracket: scan upward from 45 in steps of 5 for a sign change
    grid = [45 + 5 * i for i in range(13)]
    vals, a, b = [], None, None
    for h in grid:
        try:
            vals.append(f(h))
        except Exception as e:  # noqa: BLE001
            vals.append(None)
            print(f"  {model}: H0={h} failed: {e}")
        if len(vals) > 1 and vals[-1] is not None and vals[-2] is not None and vals[-1] * vals[-2] < 0:
            a, b = grid[len(vals) - 2], h
            break
    H0_true = brentq(f, a, b, rtol=1e-9) if a is not None else float("nan")
    th_true = theta_of(params_from(base | {"hubble": repr(H0_true)}), kind)
    r = dict(model=model, kind=kind, EFTflag=int(cp.EFTCAMB.EFTflag), H0_ini=H0_ini, th_ini=th_ini,
             rel_ini=th_ini / THETA - 1, H0_gr=H0_gr, th_hub=th_hub, H0_true=H0_true, rel_true=th_true / THETA - 1,
             dH0=H0_ini - H0_true)
    print(f"{model:24s} {kind:13s} EFTflag={r['EFTflag']} H0_ini={H0_ini:.5f} H0_EFTflag0={H0_gr:.5f} "
          f"H0_true={H0_true:.5f} | theta_ini rel={r['rel_ini']:+.3e} (hubble=H0_ini: {th_hub / THETA - 1:+.3e}) "
          f"theta(H0_true) rel={r['rel_true']:+.1e}", flush=True)
    return r


if __name__ == "__main__":
    jobs = sys.argv[1:] or ["1_EFT_GR", "4_Des_F_R_LCDM_1", "2_PEFT_wDE_wCDM_1", "5_quint_exponential",
                            "5_Kmouflage_1"]
    results = []
    with harness.slot(1):
        for job in jobs:
            if job == "fluid_w-0.8":
                results.append(check("1_EFT_GR", extra={"EFTflag": "0", "dark_energy_model": "fluid", "w": "-0.8"}))
                results[-1]["model"] = job
            elif ":" in job:
                m, k = job.split(":")
                results.append(check(m, k))
            else:
                results.append(check(job))
    tag = TREE.name + "_" + "_".join(j.replace(":", "-") for j in jobs)[:80]
    (OUT / f"indep_{tag}.json").write_text(json.dumps(results, indent=1))

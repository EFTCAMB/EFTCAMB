"""t3-gr (b): GR approach for the families not covered by b1: quintessence with an exponential potential
(EFT_lambda -> 0 flattens the potential: w -> -1), Beyond Horndeski (want_LCDM = T fixes x1 = -2 Omega_DE,
x2 = 3 Omega_DE, x3 = x4 = 0 at a_ini, the LCDM point; scan x30 = x40 -> 0 with the generic background), and the
Scaling Cubic Galileon with want_LCDM = T. Fast settings (no tensors/lensing, l_max 1500) for the model and the
EFTflag=0 reference (IntTolBoost=2, reused from b10). Reports max deviations; status/errors for failing points.
"""

import re
from concurrent.futures import ThreadPoolExecutor

from common import GR_INI, INIS, PARAMS, RUNS, fmt, harness

SUB = RUNS / "b11"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "1500", "num_redshiftwindows": "0"}
QUINT = PARAMS / "5_quint_exponential.ini"
BH = INIS / "008p5-beyond-horndeski.ini"
SCG = PARAMS / "5_scg.ini"

JOBS = {}
for lam in (1.0, 0.3, 0.1, 0.03, 0.01):
    JOBS[f"quint_exp_lam{lam:g}"] = (QUINT, {"EFT_lambda": repr(lam)})
JOBS["BH_LCDM"] = (BH, {"want_LCDM": "T"})
JOBS["BH_LCDM_x0"] = (BH, {"want_LCDM": "T", "Beyond_Horndeski_x30": "0", "Beyond_Horndeski_x40": "0"})
for s in (1.0, 1e-2, 1e-4):
    JOBS[f"BH_x34_{s:g}"] = (BH, {"Beyond_Horndeski_x30": repr(s), "Beyond_Horndeski_x40": repr(s)})
JOBS["SCG_wantLCDM"] = (SCG, {"want_LCDM": "T"})
JOBS["SCG_default"] = (SCG, {})


def one(item):
    n, (ini, ov) = item
    out = SUB / f"{n}.npz"
    if not out.with_suffix(".json").exists():
        harness.run_model("A0", ini, out, overrides={**FAST, **ov}, threads=1, timeout=540)
    return n, out


def main() -> None:
    ref = RUNS / "b10" / "REF.npz"
    if not ref.exists():
        harness.run_model("A0", GR_INI, ref, overrides=FAST, python_set={"Accuracy.IntTolBoost": 2.0}, threads=1)
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(one, JOBS.items()))
    for n, out in res.items():
        st = harness.json.loads(out.with_suffix(".json").read_text())
        log = out.with_suffix(".log").read_text(errors="replace")
        why = re.findall(r"(Model dtauda is Nan|Physical instability[^\n]*|Mathematical instability[^\n]*|"
                         r"Model has Nan[^\n]*|initial condition not found|Error in Fortran[^\n]*\n[^\n]*)", log)
        line = f"{n:20s} {st['status']:12s} {why[:1]}"
        if st["status"] == "ok":
            r = harness.compare(out, ref, lmax=1450)
            u = r["cl"].get("unlensed_scalar", {})
            line += " " + fmt(r) + f" uTT={u.get('TT', float('nan')):.1e} uEE={u.get('EE', float('nan')):.1e}"
        print(line)


if __name__ == "__main__":
    main()

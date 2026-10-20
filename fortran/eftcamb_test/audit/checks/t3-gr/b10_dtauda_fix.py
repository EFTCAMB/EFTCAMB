"""t3-gr (b): do the exact-GR points rejected with 'Model dtauda is Nan' pass the stability check on R5 (A0 +
eft_cache%grhoa2 = grhom_t*a^2 in EFTStabilityComputation, + python bridge grhoa2 fix), and do they then reproduce
LCDM? Fast settings (no tensors/lensing, l_max 1500) on both the model and the EFTflag=0 reference.
"""

import re
from concurrent.futures import ThreadPoolExecutor

from common import GR_INI, INIS, PARAMS, RUNS, WORK, fmt, harness

SUB = RUNS / "b10"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "1500", "num_redshiftwindows": "0"}
R5 = WORK / "R5"
CASES = {
    "OLalpha": (INIS / "007p6-OL-alpha.ini", {"OLmassPmodel": "0", "OLkineticitymodel": "0"}, {}),
    "Kmouflage": (PARAMS / "5_Kmouflage_1.ini", {"eps2_0": "0"}, {"omnuh2": "0.0"}),
    "Kmimic": (PARAMS / "5_Kmimic_1.ini", {"eps2_0": "0"}, {"omnuh2": "0.0"}),
    "ExtGal": (INIS / "008p7-ext-galileon.ini", {"ExtendedGalileon_S": "0", "ExtendedGalileon_q": "100"}, {}),
}
JOBS = {}
for c, (ini, gr, refov) in CASES.items():
    for t in ("A0", "R5"):
        JOBS[f"{t}_{c}"] = ((R5 if t == "R5" else "A0"), ini, {**FAST, **gr}, {})
JOBS["REF"] = ("A0", GR_INI, FAST, {"Accuracy.IntTolBoost": 2.0})
JOBS["REF_nu0"] = ("A0", GR_INI, {**FAST, "omnuh2": "0.0"}, {"Accuracy.IntTolBoost": 2.0})


def one(item):
    n, (tree, ini, ov, ps) = item
    out = SUB / f"{n}.npz"
    if not out.with_suffix(".json").exists():
        harness.run_model(tree, ini, out, overrides=ov, python_set=ps, threads=1, timeout=540)
    return n, out


with ThreadPoolExecutor(4) as ex:
    res = dict(ex.map(one, JOBS.items()))
for n, out in res.items():
    if n.startswith("REF"):
        continue
    st = harness.json.loads(out.with_suffix(".json").read_text())
    log = out.with_suffix(".log").read_text(errors="replace")
    why = re.findall(r"(Model dtauda is Nan|Physical instability[^\n]*|Mathematical instability[^\n]*|"
                     r"Model has Nan[^\n]*|Error in Fortran[^\n]*\n[^\n]*)", log)
    line = f"{n:14s} {st['status']:12s} {why[:1]}"
    if st["status"] == "ok":
        ref = res["REF_nu0" if n.split("_", 1)[1].startswith("K") else "REF"]
        line += " " + fmt(harness.compare(out, ref))
    print(line)

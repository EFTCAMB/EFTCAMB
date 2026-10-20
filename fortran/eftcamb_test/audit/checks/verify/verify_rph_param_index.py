"""Independent verification of finding rph-param-index-tensor-ode.

Oracle: EFTCAMB feedback(print_params=True). EFTCAMBRPHFeedback calls each parametrized function's own feedback in
the same block order as EFTCAMBRPHInitModelParameters / ComputeParametersNumber, and each function prints its own
(name, value) pairs with its local index, so it does not go through the model-level NMo/NKo/... offsets under test.
We compare that ordered list with EFTCAMB.param_names()/param_values() (model-level getters, 11_EFTCAMB_python.f90 ->
007p2_RPH.f90 ParameterNames/ParameterValues) for several RPH configurations, on A0 (code under audit) and on a tree
with only the NMo index changed (FIX2, which also has unrelated fixes elsewhere, none in the RPH getters).

Usage: python verify_rph_param_index.py TREE [TREE ...]
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

REPO = Path(__file__).resolve().parents[5]
PAR = REPO / "fortran/eftcamb_test/parameters"
OUT = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/work/verify/rphidx"
)

CHILD = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
e = pars.EFTCAMB
d = dict(num_params=e.num_params(), names=list(e.param_names()), values=list(e.param_values()),
         labels=list(e.param_labels()))
print("JSONSTART" + json.dumps(d) + "JSONEND", flush=True)
e.feedback(print_params=True)
"""

CASES = {
    # shipped test ini, no override: usealphaM=F, PlanckMass_ODE model 0 (0 params), Tensor_ODE model 2 (1 param)
    "KT_OmegaDE_shipped": ("3_RPH_KT_OmegaDE_1.ini", {}),
    # the finding's case: PlanckMass_ODE 1 param, Tensor_ODE 0 params
    "MK_power_massODE": ("3_RPH_MK_power_1.ini", {"RPHmassPmodel_ODE": 1, "RPHmassP_ODE0": 0.123}),
    # mislabel case: PlanckMass_ODE 1 param + Kineticity_ODE 1 param, Tensor_ODE 0 params
    "MK_power_massODE_kinODE": (
        "3_RPH_MK_power_1.ini",
        {"RPHmassPmodel_ODE": 1, "RPHmassP_ODE0": 0.123, "RPHkineticitymodel_ODE": 1, "RPHkineticity_ODE0": 0.77},
    ),
    # opposite mismatch with a trailing block: PlanckMass_ODE 0, Kineticity_ODE 1, Tensor_ODE 1
    "MK_power_kinODE_tensODE": (
        "3_RPH_MK_power_1.ini",
        {"RPHkineticitymodel_ODE": 1, "RPHkineticity_ODE0": 0.77, "RPHtensormodel_ODE": 1, "RPHtensor_ODE0": -0.05},
    ),
    # control, equal counts: PlanckMass_ODE 1, Tensor_ODE 1
    "MK_power_equal_counts": (
        "3_RPH_MK_power_1.ini",
        {"RPHmassPmodel_ODE": 1, "RPHmassP_ODE0": 0.123, "RPHtensormodel_ODE": 1, "RPHtensor_ODE0": -0.05},
    ),
    # control, usealphaM = T (alphaM_ODE branch, correct in source)
    "MK_OmegaDE_alphaM_shipped": ("3_RPH_MK_OmegaDE_1.ini", {}),
    "MB_OmegaDE_alphaM_shipped": ("3_RPH_MB_OmegaDE_1.ini", {}),
    # control, no ODE functions at all
    "MK_power_shipped": ("3_RPH_MK_power_1.ini", {}),
    "KT_power_shipped": ("3_RPH_KT_power_1.ini", {}),
    "MK_const_shipped": ("3_RPH_MK_const_1.ini", {}),
}

FB = re.compile(r"^\s*(\S+)\s*=\s*([-+0-9.EeDd]+)\s*$")


def run_case(tree: Path, ini: str, over: dict) -> dict:
    tmp = Path(tempfile.mkdtemp(prefix="vrph_", dir=OUT))
    w = tmp / "w.ini"
    lines = [f"DEFAULT({PAR / ini})", "feedback_level = 1", "EFTCAMB_write_background = F",
             f"output_root = {tmp / 'out'}"] + [f"{k} = {v}" for k, v in over.items()]
    w.write_text("\n".join(lines) + "\n")
    with harness.slot(1):
        r = subprocess.run([sys.executable, "-c", CHILD, str(tree), str(w)], capture_output=True, text=True, cwd=tmp,
                           env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=600, preexec_fn=lambda: os.nice(10))
    out = r.stdout
    m = re.search(r"JSONSTART(.*)JSONEND", out)
    if r.returncode or not m:
        return {"rc": r.returncode, "stdout": out[-2000:], "stderr": r.stderr[-2000:]}
    d = json.loads(m.group(1))
    fb = []
    in_fn = False
    raw = out[m.end():].split("Number of params", 1)[-1].splitlines()
    d["feedback_raw"] = [ln.rstrip() for ln in raw if ln.strip()]
    for ln in raw:
        if "function" in ln.lower():
            in_fn = True
            continue
        mm = FB.match(ln)
        if in_fn and mm and not ln.strip().startswith("Model") and not ln.strip().startswith("Number"):
            fb.append((mm.group(1), float(mm.group(2).replace("D", "E"))))
    d["feedback"] = fb
    d["rc"] = r.returncode
    got = list(zip(d["names"], d["values"]))
    # usealphaM=T adds one model-level parameter (RPH_M0) not printed by any function feedback: compare by name set
    fb_names = [n for n, _ in fb]
    got_f = [(n, v) for n, v in got if n in fb_names or n == ""]
    d["consistent_order"] = bool(fb) and (len(got_f) == len(fb)) and all(
        a[0] == b[0] and abs(a[1] - b[1]) < 1e-5 for a, b in zip(got_f, fb))
    d["missing_from_getters"] = [n for n in fb_names if n not in d["names"]]
    d["empty_names"] = sum(1 for n in d["names"] if n.strip() == "")
    return d


def main():
    trees = [Path(t) for t in sys.argv[1:]]
    OUT.mkdir(parents=True, exist_ok=True)
    res = {}
    for t in trees:
        tag = t.name
        res[tag] = {}
        for name, (ini, over) in CASES.items():
            d = run_case(t, ini, over)
            res[tag][name] = d
            if "names" not in d:
                print(f"[{tag}] {name}: FAILED rc={d['rc']}\n{d.get('stderr', '')}")
                continue
            print(f"[{tag}] {name}: n={d['num_params']} consistent={d['consistent_order']} "
                  f"empty={d['empty_names']} missing={d['missing_from_getters']}")
            print(f"      getters : {list(zip(d['names'], d['values']))}")
            print(f"      labels  : {d['labels']}")
            print(f"      feedback: {d['feedback']}")
            print("      raw     : " + " | ".join(d["feedback_raw"]))
    (OUT / "verify_rph_param_index.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

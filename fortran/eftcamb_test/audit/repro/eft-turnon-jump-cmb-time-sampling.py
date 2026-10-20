"""Repro (t10b-convergence): the EFT perturbation turn-on makes the potentials and the ISW source jump, and the CMB
line-of-sight time quadrature integrates that discontinuity to first order in the source time step, so the low-ell
TT of models that are not GR at EFTCAMB_turn_on_time (test inis: turn-on a = 0.01, z = 99, inside the coarse
late-time source grid dtau0 ~ 7 Mpc) depends at the percent level on TimeStepBoost / accuracy_boost and even on
unrelated redshift windows (they add source time steps).

At tau_on the scalar evolution switches from GR to the EFT equations (equations.f90 GaugeInterface_EvolveScal,
EFTturnOnTime switch, EFT_IC_type = 1 initial pi). For 2_PEFT_Omega_const_1 at k = 1e-3/Mpc Psi jumps from -0.593
to -0.452 and the ISW source EFTISW/k from +1.4e-5 to -8.1e-4 (get_scale_evolution, Part 1). The CMB sources are
sampled on State%TimeSteps (results.f90 SetTimeSteps) and integrated with trapezoid weights (cmbmain DoFlatIntegration,
TimeSteps%dpoints), so the interval that contains tau_on contributes an O(dtau0 x jump) error.

Part 2: unlensed TT at ell = 2..30, 2_PEFT_Omega_const_1 and 1_EFT_GR, redshift windows off (default), on, and
off with TimeStepBoost = 3. Expected: changes at the GR level (< 1e-3). Observed: 1-5 %.
Part 3 (if the audit tree work/t10b-convergence/knob exists): the same with AUDIT_T10B_TON=1, AUDIT_T10B_PRE=100,
which adds source time steps tightly around tau_on (dense dtau0/10 grid 100 Mpc before and 300 Mpc after):
the dependence on windows and TimeStepBoost must drop to the GR level.

Run: /Users/marco/Software/miniforge3/bin/python eft-turnon-jump-cmb-time-sampling.py   (~2 min)
"""

import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "t10b-convergence" / "repro_turnon"
KNOB = harness.AUDIT_ROOT / "work" / "t10b-convergence" / "knob"
MODEL = "2_PEFT_Omega_const_1"
ELLS = [2, 3, 4, 6, 10, 15, 20, 30]

EVOL = r'''
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_results(cp)
a = np.array([0.0099, 0.01 * (1 - 1e-6), 0.01 * (1 + 1e-6), 0.0101, 0.0105, 0.012, 0.02])
f, ev = cp.EFTCAMB.get_scale_evolution(res, [1e-3], a)
print(json.dumps({"a": a.tolist(), "Psi": ev["Psi"][0].tolist(), "ISW": (ev["EFTISW"][0] / 1e-3).tolist()}))
'''


def part1() -> None:
    tmp = Path(tempfile.mkdtemp(dir=OUT))
    ini = tmp / "w.ini"
    ini.write_text(f"DEFAULT({harness.PARAMS_DIR / (MODEL + '.ini')})\noutput_root = {tmp / 'x'}\n"
                   "EFTCAMB_write_background = F\nget_tensor_cls = F\nnum_redshiftwindows = 0\ndo_lensing = F\n"
                   "l_max_scalar = 300\nget_transfer = F\nfeedback_level = 0\n")
    (tmp / "evol.py").write_text(EVOL)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    with harness.slot(1):
        r = subprocess.run([sys.executable, str(tmp / "evol.py"), str(harness.tree_path("A0")), str(ini)], cwd=tmp,
                           env=env, capture_output=True, text=True, timeout=600)
    d = json.loads(r.stdout.strip().splitlines()[-1])
    print(f"Part 1: {MODEL}, k = 1e-3/Mpc, across the EFT turn-on a = 0.01 (get_scale_evolution, A0)")
    for a, psi, isw in zip(d["a"], d["Psi"], d["ISW"]):
        print(f"   a = {a:.7f}   Psi = {psi:+.5f}   EFTISW/k = {isw:+.4e}")


def job(args):
    tag, tree, model, over, env = args
    out = OUT / f"{model}_{tag}.npz"
    st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=over, threads=1,
                           extra_env=env, timeout=1200)
    if st["status"] != "ok":
        raise SystemExit(f"{model} {tag}: {st['status']}")
    return (model, tag), harness.load(out)["cl_unlensed_scalar"][:, 0]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    part1()
    variants = {"nowin": {"num_redshiftwindows": "0"}, "win": {},
                "ts3": {"num_redshiftwindows": "0", "TimeStepBoost": "3"}}
    jobs = [(t, "A0", m, o, {}) for m in ("1_EFT_GR", MODEL) for t, o in variants.items()]
    fix = (KNOB / "camb" / "camblib.so").exists()
    if fix:
        env = {"AUDIT_T10B_TON": "1", "AUDIT_T10B_PRE": "100"}
        jobs += [("fix_" + t, KNOB, MODEL, o, env) for t, o in variants.items()]
    with ThreadPoolExecutor(4) as ex:
        cl = dict(ex.map(job, jobs))
    print("Part 2/3: unlensed scalar TT ratio - 1 relative to the default run (windows off)")
    for m, pre in (("1_EFT_GR", ""), (MODEL, ""), (MODEL, "fix_")):
        if (m, pre + "nowin") not in cl:
            continue
        base = cl[(m, pre + "nowin")]
        for t, lab in (("win", "redshift windows on"), ("ts3", "TimeStepBoost = 3")):
            r = [cl[(m, pre + t)][l] / base[l] - 1 for l in ELLS]
            print(f"  {(pre or 'A0 ') + m:26s} {lab:20s}", " ".join(f"{l}:{x:+.2e}" for l, x in zip(ELLS, r)))
    print("expected: |ratio-1| < 1e-3 (GR level) for every row; observed: A0 2_PEFT_Omega_const_1 rows at 1e-2")


if __name__ == "__main__":
    main()

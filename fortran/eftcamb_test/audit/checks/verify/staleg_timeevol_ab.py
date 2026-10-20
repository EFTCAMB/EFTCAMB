"""Independent check of finding ic23-stale-cache-gammas on the FROZEN, uninstrumented A0 release build.

Method (different code path from the repro): CAMB's results.get_time_evolution (camb_python.f90 CAMB_TimeEvolution ->
GetOutputEvolutionFork -> initial) with 1 OpenMP thread, so modes are processed in order q(1), q(2), ... with the same
private EV.  We request q = [k, k, k] (same k three times).  Mode 1 is sacrificial (first use of EV).  The adiabatic IC
of mode 2 is computed in `initial` at taustart = 0.001/k, before EFT turns on.  If `initial` really reads Gamma3 from a
stale cache, then mode 2's delta_photon just after taustart depends on the LAST time requested for mode 1 (the last
derivs call of mode 1 fills Gamma3(a_last)).  We vary only the last requested time (a_last = 0.05, 0.3, 1) and keep all
earlier output times identical.

For 2_PEFT_gamma3_exp_1 (Omega = 0, Gamma2 = 0, Gamma3(a) = exp(0.5 a^0.1) - 1), IC type 3 gives at leading order
   delta_photon = -x^2 (2 + 3 G)/6,  x = k tau,   G = Gamma3 used by `initial`
so G_used = (6 |delta_photon| / x^2 - 2)/3 can be inverted from the output.  Correct code: G_used = Gamma3(0.01) for all
a_last.  usage: python staleg_timeevol_ab.py [model] [ic_type]
"""

import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

TREE = str(harness.TREES["A0"])
MODEL = sys.argv[1] if len(sys.argv) > 1 else "2_PEFT_gamma3_exp_1"
IC = sys.argv[2] if len(sys.argv) > 2 else "3"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "staleg"
OUT.mkdir(parents=True, exist_ok=True)

CHILD = r"""
import json, sys
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_background(pars)
k = float(sys.argv[3])
tst = 0.001 / k
early = [tst * f for f in (1.0001, 1.001, 1.01, 1.1, 2.0, 10.0)]
out = {}
for alast in (0.05, 0.3, 1.0):
    tau_last = float(res.conformal_time(1 / alast - 1))
    eta = np.array(early + [tau_last])
    ev = res.get_time_evolution([k, k, k], eta, ["delta_photon", "a"], lAccuracyBoost=1)
    out[str(alast)] = {"tau": eta.tolist(), "dg": ev[:, :, 0].tolist(), "a": ev[:, :, 1].tolist()}
print("JSON" + json.dumps(out))
"""


def run(k: float) -> dict:
    tmp = tempfile.mkdtemp(dir=OUT)
    ini = os.path.join(tmp, "w.ini")
    with open(ini, "w") as f:
        f.write(f"DEFAULT({harness.PARAMS_DIR}/{MODEL}.ini)\noutput_root = {tmp}/r\nfeedback_level = 0\n"
                f"EFTCAMB_write_background = F\nEFT_IC_type = {IC}\n")
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-c", CHILD, TREE, ini, str(k)], capture_output=True, text=True, cwd=tmp,
                           env=env, preexec_fn=lambda: os.nice(10), timeout=900)
    for line in p.stdout.splitlines():
        if line.startswith("JSON"):
            return json.loads(line[4:])
    raise RuntimeError(p.stdout[-2000:] + p.stderr[-2000:])


import math  # noqa: E402

g3 = lambda a: math.exp(0.5 * a**0.1) - 1  # noqa: E731
results = {}
for k in (0.05, 0.5):
    d = run(k)
    results[k] = d
    print(f"== {MODEL}  EFT_IC_type={IC}  k={k}  (A0 frozen release build, 1 thread, q=[k,k,k])")
    for alast, r in d.items():
        tau0 = r["tau"][0]
        x2 = (k * tau0) ** 2
        row = []
        for m in range(3):
            dg = r["dg"][m][0]
            G = (6 * abs(dg) / x2 - 2) / 3
            row.append(f"mode{m+1}: dg={dg: .6e} G_used={G: .5f}")
        print(f"  a_last={float(alast):4.2f}  " + "   ".join(row))
    print(f"  expected (correct code): G_used = Gamma3(0.01) = {g3(0.01):.5f} for every mode and every a_last")
    print(f"  stale prediction: G_used = Gamma3(a_last) = " +
          ", ".join(f"{g3(a):.5f} (a_last={a})" for a in (0.05, 0.3, 1.0)))
(OUT / f"timeevol_ab_{MODEL}_ic{IC}.json").write_text(json.dumps({str(k): v for k, v in results.items()}))

"""Repro (track t8b-integration): the theta -> H0 solve (camb/model.py set_cosmology(cosmomc_theta=...) ->
set_H0_for_theta -> scipy brentq) aborts for EFT models that are stable at the requested theta but unstable, or whose
background solver fails, at the ends of the H0 bracket.

set_H0_for_theta evaluates f(H0) = theta(H0) - theta at both bracket ends (EFTCAMB MOD lower end 1.1*100*sqrt(omega_m),
upper end 100). For EFT models, calc_background_no_thermo raises CAMBError('EFTCAMB: theory unstable' or
'EFTCAMB: background solver failed') at those trial points (SetParams runs the full EFT background + stability check
for every trial H0), and the exception propagates out of brentq. The models below are four of the 73 test models, at
their own parameters: theta computed at the ini H0 cannot be inverted back to that H0.
(Related: t8c 'theta-h0-range-eft' (lower end H0=10 with set_H0_for_theta). This shows the upper end and
set_cosmology itself, for test-suite models at their fiducial parameters.)

Expected: H0_rec == H0_ini (brentq rtol 5e-5), as obtained when the failing bracket ends are moved inwards.
Run: /Users/marco/Software/miniforge3/bin/python theta-solve-aborts-unstable-bracket-end.py
"""

import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "checks"))
try:  # hold a machine-wide CPU slot of the audit harness when available
    from harness import slot as _slot
except Exception:  # noqa: BLE001
    from contextlib import nullcontext as _slot

TREE = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
        "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
PARAMS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "parameters"))
MODELS = ["2_PEFT_Omega_power_2.ini", "3_RPH_KT_OmegaDE_1.ini", "5_quint_inverse_power.ini", "5_scg.ini"]

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
from scipy.optimize import brentq

out = {}
for ini in sys.argv[3:]:
    r = {}
    cp = camb.read_ini(sys.argv[2] + "/" + ini, no_validate=True)
    cp.EFTCAMB.EFTCAMB_feedback_level = 0
    camb.set_feedback_level(0)
    kw = dict(ombh2=cp.ombh2, omch2=cp.omch2, omk=cp.omk, mnu=None, nnu=cp.N_eff, TCMB=cp.TCMB, YHe=cp.YHe,
              omnuh2_active=cp.omnuh2, num_massive_neutrinos=max(1, int(cp.num_nu_massive)) if cp.omnuh2 > 0 else 0)
    H0_ini = cp.H0
    cp.set_cosmology(H0=H0_ini, **kw)
    theta = camb.get_background(cp, no_thermo=True).cosmomc_theta()
    r["H0_ini"], r["theta"] = H0_ini, theta
    lo = 1.1 * 100 * np.sqrt(cp.ombh2 + cp.omch2 + cp.omnuh2)

    def f(H0):
        cp.H0 = H0
        return camb.get_background(cp, no_thermo=True).cosmomc_theta() - theta

    for name, H0 in (("lower_end", lo), ("upper_end", 100.0)):
        print("probing f(H0=%.2f)" % H0, flush=True)
        try:
            r[name] = "%.2f ok" % H0 if np.isfinite(f(H0)) else "nan"
        except camb.CAMBError as e:
            r[name] = "%.2f %s" % (H0, str(e).splitlines()[-1])
    # observed: the public API
    try:
        cp.set_cosmology(cosmomc_theta=theta, **kw)
        r["observed"] = cp.H0
    except Exception as e:
        r["observed"] = repr(e).splitlines()[-1][:110]

    # expected: same brentq with failing bracket ends moved inwards (towards est_H0 = 67)
    def valid_end(a, b):
        for _ in range(40):
            try:
                return a, f(a)
            except camb.CAMBError:
                a = 0.5 * (a + b)
        raise RuntimeError("no valid bracket end")

    a, fa = valid_end(lo, 67.0)
    b, fb = valid_end(100.0, 67.0)
    r["robust_bracket"] = [a, b]
    r["expected"] = brentq(f, a, b, rtol=5e-5) if fa * fb < 0 else "no sign change"
    out[ini] = r
print("@@" + json.dumps(out))
"""


def _run(*args, **kwargs):
    with _slot(1):
        return subprocess.run(*args, preexec_fn=lambda: os.nice(10), **kwargs)


env = dict(os.environ, OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
results = {}
for model in MODELS:  # one process per model: some trial H0 values hit MpiStop (ERROR STOP) in EFTCAMB
    p = _run([sys.executable, "-c", WORKER, TREE, PARAMS, model], capture_output=True, text=True, env=env)
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@")]
    if line:
        results.update(json.loads(line[0][2:]))
    else:
        tail = [ln.strip() for ln in (p.stdout + p.stderr).splitlines() if ln.strip()][-5:]
        probe = [ln for ln in p.stdout.splitlines() if ln.startswith("probing")][-1:]
        print(f"{model}: process died (rc={p.returncode}) during {probe}; last output: {tail[-3:]}")
for ini, r in results.items():
    print(f"{ini}: H0_ini = {r['H0_ini']:.4f}, cosmomc_theta(H0_ini) = {r['theta']:.8f}")
    print(f"   f at bracket ends: lower {r['lower_end']} | upper {r['upper_end']}")
    print(f"   observed set_cosmology(cosmomc_theta) -> {r['observed']}")
    exp = r["expected"]
    exp_s = f"{exp:.4f} (rel {exp / r['H0_ini'] - 1:+.1e})" if isinstance(exp, float) else exp
    print(f"   expected (bracket ends moved inwards to {r['robust_bracket'][0]:.1f}, {r['robust_bracket'][1]:.1f})"
          f" -> H0 = {exp_s}")

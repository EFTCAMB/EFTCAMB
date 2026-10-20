"""Repro: pi'' used in the Einstein-equation factors is read from the derivative array before it is computed.

equations.f90, derivs, block "EFTCAMB MOD START: initialization of DE perturbations" (pi-field branch):

    pidotdot = ayprime(EV%w_ix + 1)          ! <- not yet written by this call
    EV%eft_cache%pidotdot = pidotdot
    call CP%EFTCAMB%model%compute_Einstein_factors(...)   ! EFTeomM and EFTeomNdot use eft_cache%pidotdot

ayprime(EV%w_ix + 1) is only filled ~270 lines later by the pi equation. In output() (CMB line-of-sight sources),
outtransf() and the python evolution fork the caller zeroes yprime first, so EFTeomM and EFTeomNdot are evaluated
with pi'' = 0; inside the Runge-Kutta stages the value is whatever the work column held (previous stage of the
previous step, or for stage 2 the error estimate of the previous step). EFTeomNdot enters EFTsigmadotdot, the ISW
source EFTISW and PsiDot; EFTeomM enters dz.

Check 1 (independent of any patch): with the A0 build (python evolution fork, same derivs code path), compare the
cached EFTeomNdot with a finite-difference derivative of the cached EFTeomN (EFTeomN does not involve pi''). The
term dropped by pi''=0 is   d EFTeomNdot = k (Gamma4 + 2 Gamma5) pi'' / (1 + Omega)   and adding it back must
restore agreement (up to the FD error, measured on pi' vs cached pi'', and the Horava Hdot inconsistency filed by
t4b, which is why the corrected residual is not at round-off).

Check 2: CMB impact. A0 vs the audit tree P built from A0 with an inert patch; AUDIT_FIX_PIDOTDOT=1 makes output()
call derivs twice so that the second call sees the correct pi'' of the first call.

Check 3 (if built): tree FX = A0 + proposed fix (recompute the Einstein factors and dz right after the pi equation in
derivs) reproduces P with AUDIT_FIX_PIDOTDOT=1 to round-off.

Run: /Users/marco/Software/miniforge3/bin/python output-pidotdot-zero.py
"""

import json
import os
import subprocess
import sys

AUDIT = (
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
    "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/"
)
TREE_A0 = AUDIT + "trees/A0"
TREE_P = AUDIT + "work/t8a-equations/P"
TREE_FX = AUDIT + "work/t8a-equations/FX"  # A0 + proposed fix (Einstein factors refreshed after pi eq.)
try:  # hold a machine-wide CPU slot (audit harness) when available; the script also runs without it
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "checks"))
    import harness  # noqa: E402

    _slot = harness.slot
except ImportError:  # pragma: no cover
    import contextlib

    _slot = lambda n=1: contextlib.nullcontext()  # noqa: E731

PARAMS = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "parameters"))

WORKER_EVOL = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.WantTensors = False
pars.InitPower.r = 0
res = camb.get_results(pars)
h0_mpc = pars.H0 / 299792.458
ks = [0.005, 0.02, 0.1]
a = np.logspace(-1.0, 0, 1500)
fields, arr = pars.EFTCAMB.get_scale_evolution(res, ks, a)
out = {}
for ik, k in enumerate(ks):
    c = arr[ik]
    tau = c["tau"]
    sl = slice(5, -5)
    fd_nd = np.gradient(c["EFTeomN"], tau, edge_order=2)
    fd_pdd = np.gradient(c["pidot"], tau, edge_order=2)
    d_nd = c["k"] * (c["EFTGamma4V"] + 2 * c["EFTGamma5V"]) * c["pidotdot"] / (1 + c["EFTOmegaV"])
    d_isw = d_nd / (c["EFTeomX"] * h0_mpc)
    norm = np.max(np.abs(fd_nd[sl]))
    out[str(k)] = {
        "pidotdot_fd_check": float(np.max(np.abs(fd_pdd[sl] - c["pidotdot"][sl])) / np.max(np.abs(fd_pdd[sl]))),
        "nd_cached_vs_fd": float(np.max(np.abs(c["EFTeomNdot"][sl] - fd_nd[sl])) / norm),
        "nd_corrected_vs_fd": float(np.max(np.abs(c["EFTeomNdot"][sl] + d_nd[sl] - fd_nd[sl])) / norm),
        "isw_missing_frac_today": float(d_isw[-1] / c["EFTISW"][-1]),
    }
print(json.dumps(out))
"""

WORKER_CL = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.WantTensors = False
pars.InitPower.r = 0
pars.DoLensing = False
res = camb.get_results(pars)
cl = res.get_unlensed_scalar_cls(CMB_unit="muK", lmax=1500)
np.save(sys.argv[3], cl)
print("ok")
"""


def run(code, args, extra_env=None):
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    env.update(extra_env or {})
    with _slot(2):
        p = subprocess.run([sys.executable, "-c", code, *args], env=env, capture_output=True, text=True,
                           cwd=PARAMS + "/..")
    if p.returncode:
        raise SystemExit(f"failed: {p.stderr[-3000:]}")
    return p.stdout.strip().splitlines()[-1]


def main():
    import numpy as np

    out_dir = AUDIT + "work/t8a-equations/repro_pidotdot"
    os.makedirs(out_dir, exist_ok=True)
    print("Check 1: cached EFTeomNdot vs finite difference of cached EFTeomN (A0, get_scale_evolution)")
    print("  expected: cached-vs-FD residual ~ FD error; observed:")
    for model in ("5_Horava_eta_4", "5_Horava_eta_3"):
        r = json.loads(run(WORKER_EVOL, [TREE_A0, os.path.join(PARAMS, model + ".ini")]))
        for k, v in r.items():
            print(f"  {model:16s} k={k:>5s}/Mpc  FD(pi')-pi''={v['pidotdot_fd_check']:.1e}  "
                  f"|Ndot_cached - FD(N)|={v['nd_cached_vs_fd']:.2e}  "
                  f"|Ndot_cached + d - FD(N)|={v['nd_corrected_vs_fd']:.2e}  "
                  f"missing/EFTISW(a=1)={v['isw_missing_frac_today']:+.2e}")
    print("Check 2: unlensed scalar Cls, A0 (pi''=0 in sources) vs P with AUDIT_FIX_PIDOTDOT=1 (correct pi'')")
    for model in ("5_Horava_eta_4", "5_Horava_eta_3"):
        ini = os.path.join(PARAMS, model + ".ini")
        fa, fp = f"{out_dir}/{model}_A0.npy", f"{out_dir}/{model}_P.npy"
        run(WORKER_CL, [TREE_A0, ini, fa])
        run(WORKER_CL, [TREE_P, ini, fp], {"AUDIT_FIX_PIDOTDOT": "1"})
        a, p = np.load(fa), np.load(fp)
        ell = np.arange(a.shape[0])
        m = ell >= 2
        dtt = (p[m, 0] - a[m, 0]) / p[m, 0]
        dte = (p[m, 3] - a[m, 3]) / np.sqrt(p[m, 0] * p[m, 1])
        i = int(np.argmax(np.abs(dtt)))
        print(f"  {model:16s} max|dTT/TT|={np.max(np.abs(dtt)):.2e} (ell={ell[m][i]})  "
              f"dTT/TT(ell=2,10,30)={dtt[0]:+.2e},{dtt[8]:+.2e},{dtt[28]:+.2e}  "
              f"max|dTE|/sqrt(TTEE)={np.max(np.abs(dte)):.2e}  max|dEE/EE|={np.max(np.abs((p[m,1]-a[m,1])/p[m,1])):.1e}"
              "  (expected 0)")
        if os.path.exists(TREE_FX + "/camb/camblib.so"):
            ff = f"{out_dir}/{model}_FX.npy"
            run(WORKER_CL, [TREE_FX, ini, ff])
            f = np.load(ff)
            print(f"  {model:16s} proposed fix FX vs P(double call): max|dTT/TT|="
                  f"{np.max(np.abs((f[m, 0] - p[m, 0]) / p[m, 0])):.1e}  (expected round-off)")


if __name__ == "__main__":
    main()

"""Repro (t10b-convergence): the pi field is integrated beyond the explicit RK45 stability limit, so the high-ell CMB
of models with fast pi oscillations depends on accuracy_boost and even on unrelated output requests.

equations.f90 GaugeInterface_ScalEv lets the scalar RK45 take steps up to 20/IntTolBoost Mpc (CAMB v2). The pi
variables are tiny in CAMB units and the mixed error norm (1/max(1,|y|)) treats them absolutely, so the error estimate
does not see a pi oscillation of frequency omega_pi^2 = (C + k^2 D1 + k^4 D2)/(A1 + k^2 A2) grow when the step exceeds
the DP45 stability limit (~3.3/omega_pi). The solution then depends on where the integrator is forced to stop
(CMB source time steps, redshift windows) and converges only slowly with the tolerance.

Check (5_Horava_SolSyst_2, test ini, redshift windows off unless stated):
  unlensed scalar TT at high ell, default accuracy vs accuracy_boost = l_accuracy_boost = l_sample_boost = 2
  (converged reference: boost 2 vs boost 3 differ by < 1e-3 for this model), and default with the two redshift
  windows of base_params.ini switched on (they only add source time steps). GR (1_EFT_GR) for comparison.
  Expected: differences at the GR level (~1e-3). Observed: ~10% at ell ~2450.
If the audit tree work/t10b-convergence/knob exists, the default run is repeated with AUDIT_T10B_PISTEP=1
(max step limited to 1.5/omega_pi/IntTolBoost while EFT is active), which must remove the discrepancy.

Run: /Users/marco/Software/miniforge3/bin/python pi-rk45-step-instability.py      (~2 min, 7 CAMB runs)
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "t10b-convergence" / "repro_pi_rk45"
KNOB = harness.AUDIT_ROOT / "work" / "t10b-convergence" / "knob"
B2 = {"accuracy_boost": "2", "l_accuracy_boost": "2", "l_sample_boost": "2"}
NOWIN = {"num_redshiftwindows": "0"}


def job(args):
    tag, tree, model, over, env = args
    out = OUT / f"{model}_{tag}.npz"
    st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=over, threads=1,
                           extra_env=env, timeout=1200)
    if st["status"] != "ok":
        raise SystemExit(f"{model} {tag}: {st['status']} {st.get('log_tail', '')[-5:]}")
    return (model, tag), harness.load(out)["cl_unlensed_scalar"][:, 0]


def main() -> None:
    jobs = []
    for m in ("1_EFT_GR", "5_Horava_SolSyst_2"):
        jobs += [("def", "A0", m, NOWIN, {}), ("b2", "A0", m, {**NOWIN, **B2}, {}), ("win", "A0", m, {}, {})]
    if (KNOB / "camb" / "camblib.so").exists():
        jobs.append(("fix", KNOB, "5_Horava_SolSyst_2", NOWIN, {"AUDIT_T10B_PISTEP": "1"}))
    with ThreadPoolExecutor(4) as ex:
        cl = dict(ex.map(job, jobs))
    ells = np.arange(1500, 2501)
    print("unlensed scalar TT, ratio - 1 (max over 1500 <= ell <= 2500, and at ell = 2000, 2450)")
    for m in ("1_EFT_GR", "5_Horava_SolSyst_2"):
        ref = cl[(m, "b2")]
        for tag, lab in (("def", "default / boost2"), ("win", "default+windows / default"), ("fix", "PISTEP fix / boost2")):
            if (m, tag) not in cl:
                continue
            base = cl[(m, "def")] if tag == "win" else ref
            r = cl[(m, tag)][ells] / base[ells] - 1
            i = int(np.argmax(np.abs(r)))
            print(f"  {m:20s} {lab:28s} max {r[i]:+.3e} at ell {ells[i]};  ell2000 {r[500]:+.3e}  ell2450 {r[950]:+.3e}")
    print("expected: all |ratio-1| ~< 1e-3 (GR level); observed: 5_Horava_SolSyst_2 default / boost2 and windows")


if __name__ == "__main__":
    main()

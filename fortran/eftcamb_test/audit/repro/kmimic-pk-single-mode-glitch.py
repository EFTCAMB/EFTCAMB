"""Repro (t10b-convergence): K-mimic test models have an isolated wrong k mode in P(k) (factor ~3) at
k ~ 0.22/Mpc, at every accuracy, and their C_ell are not converged in accuracy_boost (30-84 % at ell ~200 for
5_Kmimic_1 between default and boost 3).

k* = 0.001/tau(EFTCAMB_turn_on_time = 1e-8) ~ 0.217/Mpc is where the mode start time GetTauStart(k) = 0.001/k crosses
the EFT turn-on time: below k* the EFT initial conditions (EFT_IC_type = 1) are imposed at the mode start, above k*
after a short GR evolution at tau_on (equations.f90 GaugeInterface_EvolveScal EFT switch). For K-mimic the type-1
initial pi' is ill-conditioned (track t5-robust, repro ic1-pidot-transient-noisy-coefficients.py), so the modes next
to k* (where the two paths meet) land on a different pi branch. Which transfer-grid point is hit changes with the
k sampling, i.e. with accuracy_boost, hence the non-convergence.

Check: max |log P(k_i) - (linear interpolation in log k of log P at k_{i-1}, k_{i+1})| over 1e-3 < k < 1/Mpc for
1_EFT_GR and 5_Kmimic_1, default (windows off) and accuracy_boost = 2, on A0. Expected: ~1e-2 (BAO level, as GR).
Observed: ~1 (factor e) at k ~ 0.22. If the t5-robust tree work/t5-robust/A0_icsw exists, it is rerun with
EFT_T5_PIDOT0=1 (pi'_ini = 0, plus the EFTpiD2 typo fix), which must remove the glitch.

Run: /Users/marco/Software/miniforge3/bin/python kmimic-pk-single-mode-glitch.py   (~1 min)
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "t10b-convergence" / "repro_kmimic"
T5 = harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_icsw"
NOWIN = {"num_redshiftwindows": "0"}
B2 = {"num_redshiftwindows": "0", "accuracy_boost": "2", "l_accuracy_boost": "2", "l_sample_boost": "2"}


def glitch(d: dict) -> tuple[float, float]:
    k, p = d["pk_k"], d["pk"][0]
    m = (k > 1e-3) & (k < 1.0)
    lk, lp = np.log(k[m]), np.log(p[m])
    h1, h2 = lk[1:-1] - lk[:-2], lk[2:] - lk[1:-1]
    d2 = lp[1:-1] - (lp[:-2] * h2 + lp[2:] * h1) / (h1 + h2)
    i = int(np.argmax(np.abs(d2)))
    return float(d2[i]), float(k[m][i + 1])


def job(args):
    tag, tree, model, over, env = args
    out = OUT / f"{model}_{tag}.npz"
    st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=over, threads=1,
                           extra_env=env, timeout=1200)
    if st["status"] != "ok":
        raise SystemExit(f"{model} {tag}: {st['status']}")
    return (model, tag), harness.load(out)


def main() -> None:
    jobs = [(t, "A0", m, o, {}) for m in ("1_EFT_GR", "5_Kmimic_1") for t, o in (("def", NOWIN), ("b2", B2))]
    if (T5 / "camb" / "camblib.so").exists():
        jobs += [("t5_" + t, T5, "5_Kmimic_1", o, {"EFT_T5_PIDOT0": "1"}) for t, o in (("def", NOWIN), ("b2", B2))]
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(job, jobs))
    print("largest single-point deviation of log P(k) from its neighbours (1e-3 < k < 1/Mpc)")
    for (m, t), d in res.items():
        g, k = glitch(d)
        print(f"  {m:12s} {t:7s} {g:+.3e} at k = {k:.4f}/Mpc   sigma8 = {d['sigma8'][0]:.5f}")
    a, b = res[("5_Kmimic_1", "def")], res[("5_Kmimic_1", "b2")]
    ca, cb = a["cl_unlensed_scalar"][:, 0], b["cl_unlensed_scalar"][:, 0]
    r = ca[30:2501] / cb[30:2501] - 1
    i = int(np.argmax(np.abs(r)))
    print(f"  5_Kmimic_1 unlensed TT default/boost2 - 1: max {r[i]:+.3e} at ell {i + 30}")
    print("expected: |deviation| ~ 1e-2 for all rows (as GR); observed: A0 5_Kmimic_1 ~ 1 at k ~ 0.22")


if __name__ == "__main__":
    main()

"""Repro (t10b-convergence): the EFT_IC_type = 1 initial condition pi = -h0 E/(C + k^2 D1 + k^4 D2)
(09_EFTCAMB_IC.f90, case (1), ~l.295) diverges at the wavenumber k_c where the pi "mass + gradient" coefficient
C + k^2 D1 + k^4 D2 crosses zero at the turn-on time (C < 0: tachyonic band k < k_c). The code only guards an exact
zero. For the test model 3_RPH_MB_const_1 (turn-on a = 0.01) k_c = sqrt(-C/D1) = 1.62e-3/Mpc: the initial pi grows
like 1/(k - k_c) and changes sign, and P(k) has a pole-like spike (P/P_GR ~ 15 .. 300 next to 0.004 .. 0.08) whose
sampling depends on the transfer k grid, so P(k) and the low-ell C_ell (ell ~ k_c (tau0 - tau)) cannot converge
in accuracy_boost (default vs boost 3: P(k) differs by a factor 3 at k = 1.6e-4 .. 1.7e-3, TT by 97 % at ell 13).

Part 1: get_scale_evolution of 3_RPH_MB_const_1 just after the turn-on (a = 0.01001) for k across k_c on A0:
        C + k^2 D1, the initial pi, and delta_cdm(a=1).
Part 2: P(k)/P_GR(k) near k_c at default accuracy and accuracy_boost = 2 (windows off), A0.
        Expected: smooth in k. Observed: a spike/sign change at k ~ 1.6e-3 whose amplitude depends on the k grid.
If the audit tree work/t10b-convergence/knob exists, Part 2 is repeated with AUDIT_T10B_ICREG=0.5 (pi = pi' = 0 when
|C + k^2 D1 + k^4 D2| < 0.5 (|C| + k^2 |D1| + k^4 |D2|)) to show that the spike comes from the initial condition.

Run: /Users/marco/Software/miniforge3/bin/python ic1-pi-attractor-pole.py   (~1 min)
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

OUT = harness.AUDIT_ROOT / "work" / "t10b-convergence" / "repro_icpole"
KNOB = harness.AUDIT_ROOT / "work" / "t10b-convergence" / "knob"
MODEL = "3_RPH_MB_const_1"
KS = [0.0014, 0.0015, 0.0016, 0.00162, 0.00164, 0.0017, 0.0018, 0.0020]

EVOL = r'''
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_results(cp)
ks = json.loads(sys.argv[3])
a = np.array([0.01001, 0.1, 1.0])
f, ev = cp.EFTCAMB.get_scale_evolution(res, ks, a)
out = {k: [] for k in ("den", "C", "D1", "pi0", "clxc1")}
for i, k in enumerate(ks):
    out["C"].append(ev["EFTpiC"][i][0]); out["D1"].append(ev["EFTpiD1"][i][0])
    out["den"].append(ev["EFTpiC"][i][0] + k**2 * ev["EFTpiD1"][i][0] + k**4 * ev["EFTpiD2"][i][0])
    out["pi0"].append(ev["pi"][i][0]); out["clxc1"].append(ev["clxc"][i][-1])
print(json.dumps(out))
'''


def part1() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(dir=OUT))
    ini = tmp / "w.ini"
    ini.write_text(f"DEFAULT({harness.PARAMS_DIR / (MODEL + '.ini')})\noutput_root = {tmp / 'x'}\n"
                   "EFTCAMB_write_background = F\nget_tensor_cls = F\nnum_redshiftwindows = 0\ndo_lensing = F\n"
                   "l_max_scalar = 300\nget_transfer = F\nfeedback_level = 0\n")
    (tmp / "evol.py").write_text(EVOL)
    with harness.slot(1):
        r = subprocess.run([sys.executable, str(tmp / "evol.py"), str(harness.tree_path("A0")), str(ini), json.dumps(KS)],
                           cwd=tmp, env=dict(os.environ, OMP_NUM_THREADS="1"), capture_output=True, text=True,
                           timeout=600)
    d = json.loads(r.stdout.strip().splitlines()[-1])
    kc = np.sqrt(-d["C"][0] / d["D1"][0])
    print(f"Part 1: {MODEL} at a = 0.01001 (EFT turn-on 0.01); C = {d['C'][0]:.3e}, D1 = {d['D1'][0]:.3e}, "
          f"k_c = sqrt(-C/D1) = {kc:.4e}/Mpc")
    for k, den, pi0, c1 in zip(KS, d["den"], d["pi0"], d["clxc1"]):
        print(f"   k = {k:.5f}  C + k^2 D1 + k^4 D2 = {den:+.3e}   pi_ini = {pi0:+.3e}   delta_c(a=1) = {c1:+.3e}")


def job(args):
    tag, tree, model, over, env = args
    out = OUT / f"{model}_{tag}.npz"
    st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=over, threads=1,
                           extra_env=env, timeout=1200)
    if st["status"] != "ok":
        raise SystemExit(f"{model} {tag}: {st['status']}")
    return (model, tag), harness.load(out)


def main() -> None:
    part1()
    nowin = {"num_redshiftwindows": "0"}
    b2 = {**nowin, "accuracy_boost": "2", "l_accuracy_boost": "2", "l_sample_boost": "2"}
    jobs = [("def", "A0", m, nowin, {}) for m in ("1_EFT_GR", MODEL)] + [("b2", "A0", MODEL, b2, {})]
    if (KNOB / "camb" / "camblib.so").exists():
        jobs += [("reg_def", KNOB, MODEL, nowin, {"AUDIT_T10B_ICREG": "0.5"}),
                 ("reg_b2", KNOB, MODEL, b2, {"AUDIT_T10B_ICREG": "0.5"})]
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(job, jobs))
    g = res[("1_EFT_GR", "def")]
    print("Part 2: P(k)/P_GR(k) for 1.2e-3 < k < 2.6e-3")
    for (m, t), d in res.items():
        if m != MODEL:
            continue
        k = d["pk_k"]
        sel = (k > 1.2e-3) & (k < 2.6e-3)
        pg = np.exp(np.interp(np.log(k[sel]), np.log(g["pk_k"]), np.log(g["pk"][0])))
        print(f"  {t:8s}", " ".join(f"{kk:.5f}:{r:.3g}" for kk, r in zip(k[sel], d["pk"][0][sel] / pg)))
    print("expected: smooth ratios; observed (A0 def/b2): spike and near-zero values around k_c")


if __name__ == "__main__":
    main()

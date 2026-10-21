"""Adversarial verification of finding output-pidotdot-zero (pi'' missing in EFTeomM/EFTeomNdot).

Independent physics check: EFTLensing = sigma' + eta_k = k (Phi + Psi) and EFTISW = sigma'' + eta_k' must equal
d(EFTLensing)/dtau. EFTISW contains EFTeomNdot (with pi''), EFTLensing contains EFTeomN (no pi''), so the identity
tests whether pi'' is correctly included. Run on A0 and on FX (A0 + proposed fix), for models with Gamma4/Gamma5
(Horava, PEFT gamma4/gamma5) and controls (GR, gamma3, Omega, f(R)).

Usage: python pidotdot_verify.py evol <tree> <out.json> [models...]
"""

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import harness  # noqa: E402

PARAMS = os.path.abspath(os.path.join(HERE, "..", "..", "..", "parameters"))

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.WantTensors = False
pars.InitPower.r = 0
pars.DoLensing = False
res = camb.get_results(pars)
h0 = pars.H0 / 299792.458
import os
ks = [float(x) for x in os.environ.get("VERIFY_KS", "0.002,0.01,0.05").split(",")]
a = np.logspace(np.log10(0.05), 0, 2500)
fields, arr = pars.EFTCAMB.get_scale_evolution(res, ks, a)
out = {}
for ik, k in enumerate(ks):
    c = arr[ik]
    tau = np.array(c["tau"])
    sl = slice(10, -10)
    lens = np.array(c["EFTLensing"]); isw = np.array(c["EFTISW"])
    fd = np.gradient(lens, tau, edge_order=2)
    G45 = np.array(c["EFTGamma4V"]) + 2 * np.array(c["EFTGamma5V"])
    miss = c["k"] * G45 * np.array(c["pidotdot"]) / (1 + np.array(c["EFTOmegaV"])) / np.array(c["EFTeomX"]) / h0
    nrm = np.max(np.abs(isw[sl]))
    N = np.array(c["EFTeomN"]); Nd = np.array(c["EFTeomNdot"])
    fdN = np.gradient(N, tau, edge_order=2)
    dNd = miss * np.array(c["EFTeomX"]) * h0
    nN = max(np.max(np.abs(fdN[sl])), 1e-300)
    out[str(k)] = {
        "Ndot_vs_fdN": float(np.max(np.abs(Nd[sl] - fdN[sl])) / nN),
        "Ndot_plus_missing_vs_fdN": float(np.max(np.abs(Nd[sl] + dNd[sl] - fdN[sl])) / nN),
        "Ndot_vs_fdN_a_gt_0p3": float(np.max(np.abs((Nd - fdN)[a > 0.3][5:-10])) / nN),
        "max_abs_G4p2G5": float(np.max(np.abs(G45))),
        "cons_isw_vs_dlens": float(np.max(np.abs(isw[sl] - fd[sl])) / nrm),
        "cons_isw_plus_missing_vs_dlens": float(np.max(np.abs(isw[sl] + miss[sl] - fd[sl])) / nrm),
        "missing_over_isw_max": float(np.max(np.abs(miss[sl])) / nrm),
        "isw_a1": float(isw[-1]), "dlens_a1": float(fd[-1]), "missing_a1": float(miss[-1]),
        "pidotdot_cached_vs_fd": float(np.max(np.abs(np.gradient(np.array(c["pidot"]), tau, edge_order=2)[sl]
                                              - np.array(c["pidotdot"])[sl])) / max(np.max(np.abs(np.array(c["pidotdot"])[sl])), 1e-300)),
    }
print("RESULT " + json.dumps(out))
"""


def run_evol(tree, model):
    env = dict(os.environ, OMP_NUM_THREADS="1", PYTHONHASHSEED="0")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run(
            [sys.executable, "-c", WORKER, tree, (model if model.endswith(".ini") else os.path.join(PARAMS, model + ".ini"))],
            env=env, capture_output=True, text=True, cwd=os.path.dirname(PARAMS), timeout=1500,
            preexec_fn=lambda: os.nice(10),
        )
    for line in p.stdout.splitlines():
        if line.startswith("RESULT "):
            return json.loads(line[7:])
    return {"error": p.returncode, "tail": (p.stdout + p.stderr)[-2000:]}


def main():
    mode, tree, out = sys.argv[1], sys.argv[2], sys.argv[3]
    models = sys.argv[4:]
    from concurrent.futures import ThreadPoolExecutor

    assert mode == "evol"
    with ThreadPoolExecutor(4) as ex:
        res = dict(zip(models, ex.map(lambda m: run_evol(tree, m), models)))
    old = json.load(open(out)) if os.path.exists(out) else {}
    old.update(res)
    json.dump(old, open(out, "w"), indent=1)
    for m, r in res.items():
        if "error" in r:
            print(m, "ERROR", r["error"], r["tail"][-500:])
            continue
        for k, v in r.items():
            print(f"{m:24s} k={k:6s} |G4+2G5|max={v['max_abs_G4p2G5']:.2e} cons={v['cons_isw_vs_dlens']:.2e} "
                  f"cons+miss={v['cons_isw_plus_missing_vs_dlens']:.2e} miss/isw={v['missing_over_isw_max']:.2e} "
                  f"pdd_fd={v['pidotdot_cached_vs_fd']:.1e} Nd-fdN={v['Ndot_vs_fdN']:.2e} Nd+miss-fdN={v['Ndot_plus_missing_vs_fdN']:.2e}")


if __name__ == "__main__":
    main()

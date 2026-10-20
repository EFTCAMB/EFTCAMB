"""Independent check (verify output-pidotdot-zero): on A0 only, EFTISW should equal d(EFTLensing)/dtau.
Measure the residual, then add the analytically missing pi'' term k(G4+2G5)pi''/(1+Omega)/(X h0) and remeasure.
Also compare with an Omega-only control (G4=G5=0) where pi'' does not enter EFTeomNdot."""
import json, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import harness

TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")
INIS = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/inis/"
PARS = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/"
W = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.WantTensors = False
pars.InitPower.r = 0
res = camb.get_results(pars)
h0 = pars.H0 / 299792.458
ks = [0.005, 0.02, 0.1]
a = np.logspace(-0.7, 0, 2000)
_, arr = pars.EFTCAMB.get_scale_evolution(res, ks, a)
out = {}
sl = slice(10, -10)
for ik, k in enumerate(ks):
    c = arr[ik]; tau = c["tau"]
    fd = np.gradient(c["EFTLensing"], tau, edge_order=2)
    miss = c["k"]*(c["EFTGamma4V"] + 2*c["EFTGamma5V"])*c["pidotdot"]/(1 + c["EFTOmegaV"])/(c["EFTeomX"]*h0)
    n = np.max(np.abs(fd[sl]))
    out[str(k)] = [float(np.max(np.abs(c["EFTISW"][sl]-fd[sl]))/n), float(np.max(np.abs(c["EFTISW"][sl]+miss[sl]-fd[sl]))/n),
                   float(np.max(np.abs(miss[sl]))/n)]
print(json.dumps(out))
"""
for ini in [INIS + "verify_v2/verify_pidd_omega_g5m.ini", INIS + "verify_v2/verify_pidd_omega_g4m.ini", INIS + "verify_v2/verify_pidd_omega_only.ini",
            PARS + "1_EFT_GR.ini"]:
    env = dict(os.environ, OMP_NUM_THREADS="1"); env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-c", W, TREE, ini], env=env, capture_output=True, text=True,
                           cwd="/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/pidotdot_v2")
    if p.returncode:
        print(os.path.basename(ini), "FAIL", p.stderr[-800:]); continue
    r = json.loads(p.stdout.strip().splitlines()[-1])
    for k, (r0, r1, m) in r.items():
        print(f"{os.path.basename(ini):28s} k={k:6s} |ISW-dLens/dtau|={r0:.2e}  with missing pi'' term={r1:.2e}  |missing|={m:.2e}")

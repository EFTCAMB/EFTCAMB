"""Repro (track t8b-integration): changing the dark-energy model of a CAMBparams that holds an active EFTCAMB model
(EFTflag != 0) is silently ignored.

camb/model.py set_classes (EFTCAMB MOD) says "When dark_energy_model \\= EFTCAMB EFTCAMB is initialized with
EFTflag = 0", but the elif branch only replaces CAMBparams.DarkEnergy and leaves CAMBparams.EFTCAMB untouched;
CAMBparams.set_dark_energy does not touch EFTCAMB either. Since dtauda/perturbations use the EFT model whenever
EFTflag != 0, the newly selected CAMB dark energy has no effect, with no warning.

Model: 2_PEFT_wDE_wCDM_1 (pure EFT, wCDM background w = -0.9). Then switch to CAMB's PPF fluid with w = -0.7.
Expected: H(z=1) equals that of a plain CAMB w = -0.7 model (EFTflag 0). Observed: still the EFT model's H(z=1).
Run: /Users/marco/Software/miniforge3/bin/python set-dark-energy-ignored-with-active-eft.py
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

TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
        "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
PARAMS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "parameters"))

WORKER = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb

def H1(cp):
    return float(camb.get_background(cp, no_thermo=True).hubble_parameter(1.0))

out = {}
cp = camb.read_ini(sys.argv[2] + "/2_PEFT_wDE_wCDM_1.ini", no_validate=True)
cp.EFTCAMB.EFTCAMB_feedback_level = 0
out["EFT model"] = (int(cp.EFTCAMB.EFTflag), H1(cp))
cp.set_dark_energy(w=-0.7, dark_energy_model="ppf")
out["after set_dark_energy(w=-0.7, 'ppf')"] = (int(cp.EFTCAMB.EFTflag), H1(cp))
cp = camb.read_ini(sys.argv[2] + "/2_PEFT_wDE_wCDM_1.ini", no_validate=True)
cp.EFTCAMB.EFTCAMB_feedback_level = 0
cp.set_classes(dark_energy_model="ppf")
cp.DarkEnergy.set_params(w=-0.7)
out["after set_classes('ppf') + DarkEnergy.set_params(w=-0.7)"] = (int(cp.EFTCAMB.EFTflag), H1(cp))
ref = camb.set_params(H0=cp.H0, ombh2=cp.ombh2, omch2=cp.omch2, omk=cp.omk, mnu=0.0, num_massive_neutrinos=0,
                      nnu=cp.N_eff, TCMB=cp.TCMB, YHe=cp.YHe, w=-0.7, dark_energy_model="ppf")
out["plain CAMB ppf w=-0.7 (expected)"] = (int(ref.EFTCAMB.EFTflag), H1(ref))
print("@@" + json.dumps(out))
"""

env = dict(os.environ, OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
with _slot(1):
    p = subprocess.run([sys.executable, "-c", WORKER, TREE, PARAMS], capture_output=True, text=True, env=env,
                       preexec_fn=lambda: os.nice(10))
line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@")]
if not line:
    print(p.stdout[-3000:], p.stderr[-3000:])
    sys.exit(1)
res = json.loads(line[0][2:])
expected = res["plain CAMB ppf w=-0.7 (expected)"][1]
for k, (flag, h) in res.items():
    print(f"   {k:58s} EFTflag={flag}  H(z=1) = {h:.4f} km/s/Mpc  (vs expected {h / expected - 1:+.2e})")

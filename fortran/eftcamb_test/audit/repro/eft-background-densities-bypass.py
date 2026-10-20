"""Repro (track t8b-integration): CAMBdata.get_background_densities / get_Omega('de') / get_dark_energy_rho_w
ignore the EFT background when EFTflag != 0 (Fortran CAMBdata_GetBackgroundDensities and
CAMBdata_DarkEnergyStressEnergy call CP%DarkEnergy%BackgroundDensityAndPressure, i.e. the never-initialised
w=-1 placeholder fluid), so they return LCDM values inconsistent with H(z) from the same results object.

Observed: H reconstructed from the returned densities vs results.hubble_parameter (dtauda -> EFT).
Expected: agreement to round-off (as for EFTflag=0).
Run: /Users/marco/Software/miniforge3/bin/python eft-background-densities-bypass.py
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


def _run(*args, **kwargs):
    with _slot(1):
        return subprocess.run(*args, preexec_fn=lambda: os.nice(10), **kwargs)


TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
        "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
PARAMS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "parameters")

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb

c = 299792.458
out = {}
for ini in sys.argv[3:]:
    cp = camb.read_ini(sys.argv[2] + "/" + ini, no_validate=True)
    res = camb.get_background(cp)
    z = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
    a = 1 / (1 + z)
    H = np.asarray(res.hubble_parameter(z))
    d = res.get_background_densities(a, vars=["tot", "de"])
    H_dens = np.sqrt(d["tot"] / 3) / a**2 * c
    om_de_api = np.asarray(res.get_Omega("de", z))
    om_de_true = 1 - (d["tot"] - d["de"]) / (3 * (H / c) ** 2 * a**4)
    rho, w = res.get_dark_energy_rho_w(a)
    out[ini] = dict(EFTflag=int(cp.EFTCAMB.EFTflag), z=z.tolist(), H=H.tolist(), H_dens=H_dens.tolist(),
                    om_api=om_de_api.tolist(), om_true=om_de_true.tolist(), w_api=np.asarray(w).tolist())
print("@@" + json.dumps(out))
"""

inis = ["1_EFT_GR.ini", "2_PEFT_wDE_wCDM_1.ini", "4_Des_mc5e_CPL_1.ini", "5_Kmouflage_1.ini"]
env = dict(os.environ, OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
p = _run([sys.executable, "-c", WORKER, TREE, os.path.abspath(PARAMS)] + inis,
                   capture_output=True, text=True, env=env)
line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@")]
if not line:
    print(p.stdout[-2000:], p.stderr[-2000:])
    sys.exit(1)
res = json.loads(line[0][2:])
for ini, r in res.items():
    rel = max(abs(x / y - 1) for x, y in zip(r["H_dens"], r["H"]))
    print(f"{ini}  (EFTflag={r['EFTflag']})")
    print("   z                         ", r["z"])
    print("   H(z) hubble_parameter      ", ["%.3f" % x for x in r["H"]])
    print("   H(z) from get_background_densities (observed)", ["%.3f" % x for x in r["H_dens"]])
    print("   Omega_de get_Omega (observed)  ", ["%.4f" % x for x in r["om_api"]])
    print("   Omega_de implied by H (expected)", ["%.4f" % x for x in r["om_true"]])
    print("   w_de get_dark_energy_rho_w      ", ["%.3f" % x for x in r["w_api"]])
    print(f"   max |H_dens/H - 1| = {rel:.3e}  (expected ~1e-15)")

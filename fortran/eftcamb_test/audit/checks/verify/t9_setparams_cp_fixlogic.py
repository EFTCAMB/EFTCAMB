"""Runtime check of the proposed fix logic for set-params-cp-drops-eft-model (no rebuild, no repo edit).

EFTCAMB.initialize_parameters is monkeypatched to merge the parameters already read by this EFTCAMB object (when
EFTflag != 0) under this call's keywords, which is what the proposed EFT_do_set change does (the diff uses
read_parameters(raw=True); here the autoconverted values are re-stringified, enough for these models).
Cases, starting from camb.read_ini(2_PEFT_Omega_const_1.ini) as in the repro:
  A: set_params(cp=p, H0=p.H0) vs control read_ini + set_cosmology(H0)
  B: set_params(cp=p, EFTOmega0=0.02) vs fresh read_ini with EFTOmega0 = 0.02 + set_cosmology(H0)? (no: B only
     touches EFT, so control = read_ini with EFTOmega0=0.02 override)
  D: set_params(cp=p, EFTflag=0) -> must switch EFT off (GR), as an explicit request

Run: /Users/marco/Software/miniforge3/bin/python t9_setparams_cp_fixlogic.py
"""

import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import harness  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, "../../../../.."))
INI = os.path.join(REPO, "fortran/eftcamb_test/parameters/2_PEFT_Omega_const_1.ini")
TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
        "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
OUT = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
       "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/t9_setparams_cp")

CHILD = r"""
import sys, os, json
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
from camb import eftcamb
ini, tmp = sys.argv[2], sys.argv[3]
def wrap(extra=''):
    f = os.path.join(tmp, 'model%d.ini' % len(extra))
    open(f, 'w').write(f'DEFAULT({ini})\nl_max_scalar = 1200\nget_tensor_cls = F\nget_transfer = F\n'
                       f'EFTCAMB_write_background = F\nfeedback_level = 0\n{extra}')
    return f
def tt(p):
    return camb.get_results(p).get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)['total'][2:1200, 0]
orig = eftcamb.EFTCAMB.initialize_parameters
def patched(self, cp, params, print_header=True):
    merged = dict(cp.EFTCAMB.read_parameters()) if cp.EFTCAMB.EFTflag != 0 else {}
    merged.update(params)
    return orig(self, cp, merged, print_header=print_header)
eftcamb.EFTCAMB.initialize_parameters = patched
out = {}
q = camb.read_ini(wrap()); q.set_cosmology(H0=q.H0); tc = tt(q)
p = camb.read_ini(wrap()); camb.set_params(cp=p, H0=p.H0)
out['A'] = dict(EFTflag=int(p.EFTCAMB.EFTflag), maxrel=float(np.max(np.abs(tt(p) / tc - 1))))
q = camb.read_ini(wrap('EFTOmega0 = 0.02\n')); tc = tt(q)
p = camb.read_ini(wrap()); camb.set_params(cp=p, EFTOmega0=0.02)
out['B'] = dict(EFTflag=int(p.EFTCAMB.EFTflag), values=list(p.EFTCAMB.param_values()),
                maxrel=float(np.max(np.abs(tt(p) / tc - 1))))
p = camb.read_ini(wrap()); camb.set_params(cp=p, EFTflag=0)
out['D'] = dict(EFTflag=int(p.EFTCAMB.EFTflag))
print('JSON' + json.dumps(out))
"""


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with tempfile.TemporaryDirectory(dir=OUT) as tmp, harness.slot(1):
        p = subprocess.run([sys.executable, "-c", CHILD, TREE, INI, tmp], capture_output=True, text=True, env=env,
                           cwd=tmp, preexec_fn=lambda: os.nice(10), timeout=1200)
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("JSON")]
    if not line:
        raise SystemExit(f"failed rc={p.returncode}:\n{p.stdout[-3000:]}\n{p.stderr[-3000:]}")
    res = json.loads(line[0][4:])
    print(json.dumps(res, indent=1))
    with open(os.path.join(OUT, "t9_setparams_cp_fixlogic.json"), "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()

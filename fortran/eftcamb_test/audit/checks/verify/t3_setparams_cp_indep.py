"""Independent check of finding t3-set-params-cp-resets-eft-model.

Different path from the repro: the EFT model comes from camb.read_ini(<test ini>) (Fortran ini reader), not from
set_params keywords, and the reuse call passes either no keywords at all or only a plain attribute (WantTensors),
so no CAMB setter resets anything. In upstream CAMB both calls are no-ops on cp. Checks per model:
  ref      : sigma8 / TT of the pars as read
  empty    : set_params(cp=pars)                    -> EFTflag, model, sigma8
  attr     : set_params(cp=pars, WantTensors=False) -> EFTflag, model, sigma8
  B0_GR    : plain CAMB (B0) read_ini of the same ini (ignores EFT keys)  = independent GR reference
  fix      : set_params(cp=pars) with EFTCAMB.initialize_parameters monkeypatched to skip when the dict has no EFT
             key (the proposed fix, applied at runtime) -> should equal ref bit for bit
Control on B0: set_params(cp=pars) and set_params(cp=pars, WantTensors=False) leave sigma8 bit-identical.

Run: /Users/marco/Software/miniforge3/bin/python t3_setparams_cp_indep.py
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
PARS = os.path.join(REPO, "fortran/eftcamb_test/parameters")
TREES = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/")
OUT = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
       "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/t3_setparams_cp")
MODELS = ["4_Des_F_R_LCDM_1", "5_quint_exponential", "2_PEFT_gamma1_const_1"]

CHILD = r"""
import sys, json, os
tree, tag, ini = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, tree)
import numpy as np
import camb

def prep(p):
    p.WantTensors = False
    p.InitPower.r = 0.0
    p.set_for_lmax(1500, lens_potential_accuracy=0)
    p.set_matter_power(redshifts=[0.0], kmax=2.0)
    p.WantTransfer = True
    return p

def obs(p):
    r = camb.get_results(p)
    tt = r.get_cmb_power_spectra(CMB_unit='muK', raw_cl=True, lmax=1500)['total'][2:1501, 0]
    return float(r.get_sigma8_0()), tt

def state(p):
    if tag == 'B0':
        return {}
    return {'EFTflag': int(p.EFTCAMB.EFTflag), 'model': p.EFTCAMB.model_name().strip()}

out = {}
def read():
    return prep(camb.read_ini(ini))

p = read(); s, tt0 = obs(p); out['ref'] = dict(sigma8=s, **state(p))
p = read(); p = camb.set_params(cp=p); s, tt = obs(p)
out['empty'] = dict(sigma8=s, tt_maxrel=float(np.max(np.abs(tt / tt0 - 1))), **state(p))
p = read(); p = camb.set_params(cp=p, WantTensors=False); s, tt = obs(p)
out['attr'] = dict(sigma8=s, tt_maxrel=float(np.max(np.abs(tt / tt0 - 1))), **state(p))
if tag == 'A0':
    from camb import eftcamb
    orig = eftcamb.EFTCAMB.initialize_parameters
    calls = []
    def patched(self, cp, params, print_header=True):
        keys = set(self.read_parameters()) | {'EFTflag'}
        hit = keys & set(params)
        calls.append(sorted(hit))
        if hit:
            return orig(self, cp, params, print_header=print_header)
    eftcamb.EFTCAMB.initialize_parameters = patched
    p = read(); p = camb.set_params(cp=p); s, tt = obs(p)
    out['fix'] = dict(sigma8=s, tt_maxrel=float(np.max(np.abs(tt / tt0 - 1))), calls=calls, **state(p))
    eftcamb.EFTCAMB.initialize_parameters = orig
print('JSON' + json.dumps(out))
"""


def run(tree: str, ini: str) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with tempfile.TemporaryDirectory(dir=OUT) as tmp, harness.slot(1):
        wrap = os.path.join(tmp, "model.ini")
        with open(wrap, "w") as f:
            f.write(f"DEFAULT({ini})\nEFTCAMB_write_background = F\nfeedback_level = 0\n")
        p = subprocess.run([sys.executable, "-c", CHILD, TREES + tree, tree, wrap], capture_output=True, text=True,
                           env=env, cwd=tmp, preexec_fn=lambda: os.nice(10), timeout=1200)
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("JSON")]
    if not line:
        raise SystemExit(f"{tree} {ini} failed:\n{p.stdout[-3000:]}\n{p.stderr[-3000:]}")
    return json.loads(line[0][4:])


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    res = {}
    for m in MODELS:
        ini = os.path.join(PARS, m + ".ini")
        a0, b0 = run("A0", ini), run("B0", ini)
        res[m] = {"A0": a0, "B0": b0}
        r, e, at, fx = a0["ref"], a0["empty"], a0["attr"], a0["fix"]
        print(f"== {m}")
        print(f"  A0 ref (read_ini)              : EFTflag={r['EFTflag']} model={r['model']!r} sigma8={r['sigma8']:.6f}")
        print(f"  A0 set_params(cp=p)            : EFTflag={e['EFTflag']} model={e['model']!r} sigma8={e['sigma8']:.6f}"
              f" max|dTT/TT|={e['tt_maxrel']:.3e}")
        print(f"  A0 set_params(cp=p,WantTensors): EFTflag={at['EFTflag']} sigma8={at['sigma8']:.6f}"
              f" max|dTT/TT|={at['tt_maxrel']:.3e}")
        print(f"  B0 read_ini (GR reference)     : sigma8={b0['ref']['sigma8']:.6f}"
              f"  rel diff vs A0 empty = {e['sigma8'] / b0['ref']['sigma8'] - 1:.2e}")
        print(f"  B0 set_params(cp=p) control    : sigma8={b0['empty']['sigma8']:.6f} max|dTT/TT|="
              f"{b0['empty']['tt_maxrel']:.3e}; attr: max|dTT/TT|={b0['attr']['tt_maxrel']:.3e}")
        print(f"  A0 runtime fix (skip w/o keys) : EFTflag={fx['EFTflag']} sigma8={fx['sigma8']:.6f}"
              f" max|dTT/TT|={fx['tt_maxrel']:.3e} calls={fx['calls']}")
    with open(os.path.join(OUT, "t3_setparams_cp_indep.json"), "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()

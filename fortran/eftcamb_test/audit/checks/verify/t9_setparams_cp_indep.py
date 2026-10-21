"""Independent check of finding set-params-cp-drops-eft-model (t9-wrapper).

Different path from the repro: the EFT model is created by camb.set_params(**keywords) (no read_ini), and the reuse
call changes nothing in the CAMB part at all (set_params(cp=p, lmax=<same lmax>) only calls set_for_lmax with the
value p already has), so the expected result is bit-identical to p. Two models (Designer f(R), Pure EFT).
  ref      : the pars as created
  A_same   : set_params(cp=p, lmax=same)                -> EFTflag, model, sigma8, max|dTT/TT| vs ref
  GR       : set_params(<same CAMB keywords>) without EFT keys -> compare A_same to it (is A_same exactly GR?)
  B_one    : set_params(cp=p, <one EFT parameter changed>) -> exception? state afterwards
  B_direct : p.EFTCAMB.initialize_parameters(p, {<one EFT key>}) -> EFTflag and read keys (Fortran-level reset)
  C_full   : set_params(cp=p, <all EFT keys, one changed>) -> works (the workaround), and equals a fresh build
Upstream control (B0): p with DarkEnergyPPF w=-0.9; set_params(cp=p, lmax=same) keeps w and TT bit-identical.
Pre-rebase (B1): A_same and B_one, to see whether the behaviour predates the rebase.

Run: /Users/marco/Software/miniforge3/bin/python t9_setparams_cp_indep.py
"""

import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import harness  # noqa: E402

TREES = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/")
OUT = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
       "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/t9_setparams_cp")

CHILD = r"""
import sys, json
tree, tag = sys.argv[1], sys.argv[2]
sys.path.insert(0, tree)
import numpy as np
import camb

COSMO = dict(H0=67.5, ombh2=0.0224, omch2=0.12, YHe=0.245, As=2.1e-9, ns=0.965, lmax=1200, WantTransfer=True,
             kmax=2.0)
MODELS = {
    'DesFR': (dict(EFTflag=3, DesignerEFTmodel=1, EFTwDE=0, EFTB0=0.1), 'EFTB0', 0.2),
    'PEFT': (dict(EFTflag=1, PureEFTmodel=1, PureEFTHorndeski=False, EFTwDE=0, PureEFTmodelOmega=1,
                  PureEFTmodelGamma1=0, PureEFTmodelGamma2=0, PureEFTmodelGamma3=0, PureEFTmodelGamma4=0,
                  PureEFTmodelGamma5=0, PureEFTmodelGamma6=0, EFTOmega0=0.1), 'EFTOmega0', 0.02),
}

def obs(p):
    r = camb.get_results(p)
    tt = r.get_cmb_power_spectra(CMB_unit='muK', raw_cl=True, lmax=1200)['total'][2:1201, 0]
    return float(r.get_sigma8_0()), tt

def state(p):
    if tag == 'B0':
        return {}
    d = {'EFTflag': int(p.EFTCAMB.EFTflag), 'model': p.EFTCAMB.model_name().strip()}
    return d

def rel(a, b):
    return float(np.max(np.abs(a / b - 1)))

out = {}
if tag == 'B0':
    p = camb.set_params(dark_energy_model='DarkEnergyPPF', w=-0.9, **COSMO)
    s0, tt0 = obs(p)
    p = camb.set_params(cp=p, lmax=1200)
    s1, tt1 = obs(p)
    out['B0_control'] = dict(w_after=float(p.DarkEnergy.w), de=type(p.DarkEnergy).__name__, sigma8_ref=s0,
                             sigma8_after=s1, tt_maxrel=rel(tt1, tt0))
else:
    g_s, g_tt = obs(camb.set_params(**COSMO))
    out['GR'] = dict(sigma8=g_s)
    for name, (eft, key, newval) in MODELS.items():
        o = {}
        p = camb.set_params(**COSMO, **eft)
        s0, tt0 = obs(p)
        o['ref'] = dict(sigma8=s0, tt_maxrel_vs_GR=rel(tt0, g_tt), **state(p))
        camb.set_params(cp=p, lmax=1200)
        s1, tt1 = obs(p)
        o['A_same'] = dict(sigma8=s1, tt_maxrel_vs_ref=rel(tt1, tt0), tt_maxrel_vs_GR=rel(tt1, g_tt),
                           sigma8_minus_GR=s1 - g_s, **state(p))
        p = camb.set_params(**COSMO, **eft)
        try:
            camb.set_params(cp=p, **{key: newval})
            exc = 'none'
        except Exception as e:
            exc = f'{type(e).__name__}: {e}'
        o['B_one'] = dict(exception=exc, **state(p))
        if tag == 'A0':
            p = camb.set_params(**COSMO, **eft)
            p.EFTCAMB.initialize_parameters(p, {key: newval}, print_header=False)
            o['B_direct'] = dict(read_keys=sorted(p.EFTCAMB.read_parameters()), **state(p))
            eft2 = dict(eft, **{key: newval})
            p = camb.set_params(**COSMO, **eft)
            camb.set_params(cp=p, **eft2)
            s2, tt2 = obs(p)
            sf, ttf = obs(camb.set_params(**COSMO, **eft2))
            o['C_full'] = dict(values=list(p.EFTCAMB.param_values()), sigma8=s2, sigma8_fresh=sf,
                               tt_maxrel_vs_fresh=rel(tt2, ttf), tt_maxrel_vs_ref=rel(tt2, tt0), **state(p))
        out[name] = o
print('JSON' + json.dumps(out))
"""


def run(tree: str) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with tempfile.TemporaryDirectory(dir=OUT) as tmp, harness.slot(1):
        p = subprocess.run([sys.executable, "-c", CHILD, TREES + tree, tree], capture_output=True, text=True,
                           env=env, cwd=tmp, preexec_fn=lambda: os.nice(10), timeout=1200)
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("JSON")]
    if not line:
        raise SystemExit(f"{tree} failed rc={p.returncode}:\n{p.stdout[-3000:]}\n{p.stderr[-3000:]}")
    return json.loads(line[0][4:])


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    trees = sys.argv[1:] or ["A0", "B0", "B1"]
    res = {}
    for t in trees:
        res[t] = run(t)
        print(f"===== {t}")
        print(json.dumps(res[t], indent=1))
    with open(os.path.join(OUT, "t9_setparams_cp_indep.json"), "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()

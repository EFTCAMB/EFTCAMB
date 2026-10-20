"""Repro (t3-gr): camb.set_params(cp=existing_params, ...) silently switches an EFT model back to GR (EFTflag=0).

camb/camb.py set_params (EFTCAMB MOD "dedicated setter for the initialization of EFTCAMB") calls
cp.EFTCAMB.initialize_parameters(cp, params) unconditionally, unlike CAMB's do_set, which calls a setter only when
one of its arguments is in params. initialize_parameters re-initialises EFTCAMB from the params dict, so when the
dict has no EFT keys (e.g. updating only H0 on an existing CAMBparams, a documented use of cp=) EFTflag becomes 0
and the results are GR, with no warning. The same pattern for the CAMB dark-energy model keeps w (control rows).

Expected: EFTflag and sigma8 of the f(R) model unchanged after set_params(cp=pars, H0=68) (as w is for the fluid).
Run: /Users/marco/Software/miniforge3/bin/python t3-set-params-cp-resets-eft-model.py
"""

import json
import os
import subprocess
import sys

TREES = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/")

WORKER = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
out = {}
base = dict(ombh2=0.0224, omch2=0.12, mnu=0.06, As=2.1e-9, ns=0.965, tau=0.055, WantCls=False, WantTransfer=True,
            kmax=1.0)

def s8(pars):
    return float(camb.get_results(pars).get_sigma8_0())

# control: CAMB fluid dark energy is kept when cp= is reused with only H0
pars = camb.set_params(H0=67.5, w=-0.9, dark_energy_model="fluid", **base)
pars2 = camb.set_params(cp=pars, H0=68.0)
out["fluid_w_after"] = float(pars2.DarkEnergy.w)
if sys.argv[2] == "A0":
    eft = {"EFTflag": 3, "DesignerEFTmodel": 1, "EFTwDE": 0, "EFTB0": 0.1}
    pars = camb.set_params(H0=68.0, **base, **eft)
    out["EFTflag_new"] = int(pars.EFTCAMB.EFTflag)
    out["sigma8_fR_H068"] = s8(pars)
    pars = camb.set_params(H0=67.5, **base, **eft)
    pars = camb.set_params(cp=pars, H0=68.0)
    out["EFTflag_after_cp_reuse"] = int(pars.EFTCAMB.EFTflag)
    out["sigma8_after_cp_reuse"] = s8(pars)
    # GR through the same cp-reuse path (set_cosmology(H0=68) also resets ombh2 etc. to their defaults)
    pars = camb.set_params(H0=67.5, **base)
    out["sigma8_GR_after_cp_reuse"] = s8(camb.set_params(cp=pars, H0=68.0))
print("JSON" + json.dumps(out))
"""


def run(tree: str) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", WORKER, TREES + tree, tree], capture_output=True, text=True, env=env)
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("JSON")]
    if not line:
        raise SystemExit(f"{tree} failed:\n{p.stdout[-2000:]}\n{p.stderr[-2000:]}")
    return json.loads(line[0][4:])


def main() -> None:
    b0, a0 = run("B0"), run("A0")
    print(f"control (B0, CAMB v2): fluid w after set_params(cp=pars, H0=68) = {b0['fluid_w_after']}  (set -0.9)")
    print(f"control (A0)         : fluid w after set_params(cp=pars, H0=68) = {a0['fluid_w_after']}  (set -0.9)")
    print(f"A0 f(R) B0=0.1, new params at H0=68     : EFTflag={a0['EFTflag_new']}  sigma8={a0['sigma8_fR_H068']:.6f}")
    print(f"A0 f(R) B0=0.1, set_params(cp=, H0=68)  : EFTflag={a0['EFTflag_after_cp_reuse']}  "
          f"sigma8={a0['sigma8_after_cp_reuse']:.6f}")
    print(f"A0 GR (EFTflag=0), set_params(cp=, H0=68): sigma8={a0['sigma8_GR_after_cp_reuse']:.6f}")
    print("expected: EFTflag=3 and sigma8 equal to the 'new params' row after cp reuse; observed: GR values")


if __name__ == "__main__":
    main()

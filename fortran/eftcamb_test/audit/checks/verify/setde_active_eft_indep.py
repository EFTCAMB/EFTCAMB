"""Independent verification of finding set-dark-energy-ignored-with-active-eft.

Claim: on a CAMBparams with an active EFT model (EFTflag != 0), set_classes(dark_energy_model=...),
set_dark_energy(...) and set_dark_energy_w_a(...) change only CAMBparams.DarkEnergy; EFTflag stays on and the
Fortran (dtauda, perturbations) ignores CP%DarkEnergy, so the change has no effect and no warning is raised.

Independent method: instead of comparing with a separately constructed plain-CAMB model (as the repro does, with
different neutrino settings), the "expected" value is computed on the *same* params object with EFTflag forced to 0
(cp.EFTCAMB.EFTflag = 0) before the dark-energy call, so every other parameter is identical. Checks both the
background (H(z=1)) and the full CMB TT spectrum (get_results), plus the camb.set_params and set_params_cosmomc
paths, the default EFTflag of a fresh CAMBparams, and the pre-rebase tree B1.

Run: /Users/marco/Software/miniforge3/bin/python setde_active_eft_indep.py
"""

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import harness  # noqa: E402

TREES = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/")
INI = os.path.abspath(os.path.join(HERE, "..", "..", "..", "parameters", "2_PEFT_wDE_wCDM_1.ini"))

WORKER = r"""
import sys, json, warnings
import numpy as np
sys.path.insert(0, sys.argv[1])
sys.path.insert(0, sys.argv[3])
import camb, harness
INI = sys.argv[2]
ini_keys = harness.parse_ini(INI)
warnings.simplefilter("always")

def load():
    cp = camb.read_ini(INI, no_validate=True)
    if cp.EFTCAMB.EFTflag == 0 and int(ini_keys.get("EFTflag", "0")) != 0:   # B1 fallback
        cp.EFTCAMB.initialize_parameters(cp, ini_keys, print_header=False)
    cp.EFTCAMB.EFTCAMB_feedback_level = 0
    cp.set_for_lmax(1500, lens_potential_accuracy=0)
    return cp

def H1(cp):
    return float(camb.get_background(cp, no_thermo=True).hubble_parameter(1.0))

def TT(cp):
    r = camb.get_results(cp)
    return r.get_cmb_power_spectra(cp, CMB_unit="muK", raw_cl=True)["total"][2:1501, 0]

out = {"tree": sys.argv[1].rstrip("/").split("/")[-1]}
fresh = camb.CAMBparams()
out["fresh CAMBparams EFTflag"] = int(fresh.EFTCAMB.EFTflag) if fresh.EFTCAMB is not None else None

cp = load()
out["EFTflag after read_ini"] = int(cp.EFTCAMB.EFTflag)
out["DE after read_ini"] = [type(cp.DarkEnergy).__name__, float(cp.DarkEnergy.w)]
h_eft = H1(cp); tt_eft = TT(cp)
out["H1 EFT model"] = h_eft

# sanity: EFT wCDM model (EFTw0=-0.8) equals CAMB fluid w=-0.8 on identical params with EFTflag=0
c = load(); c.EFTCAMB.EFTflag = 0; c.set_dark_energy(w=-0.8, dark_energy_model="fluid")
out["H1 EFTflag=0 fluid w=-0.8 (sanity, should equal EFT)"] = H1(c)

# expected after switching to ppf w=-0.7: same params, EFTflag forced 0
c = load(); c.EFTCAMB.EFTflag = 0; c.set_dark_energy(w=-0.7, dark_energy_model="ppf")
h_exp = H1(c); tt_exp = TT(c)
out["H1 expected (EFTflag=0, ppf w=-0.7)"] = h_exp

cases = {}
def record(name, c, w_caught, do_tt=False):
    d = {"EFTflag": int(c.EFTCAMB.EFTflag), "DE": [type(c.DarkEnergy).__name__, float(c.DarkEnergy.w)],
         "H1": H1(c), "warnings": [str(x.message) for x in w_caught]}
    d["H1 rel vs EFT"] = d["H1"] / h_eft - 1
    d["H1 rel vs expected"] = d["H1"] / h_exp - 1
    if do_tt:
        tt = TT(c)
        d["TT bit-identical to EFT"] = bool(np.array_equal(tt, tt_eft))
        d["TT max rel vs EFT"] = float(np.max(np.abs(tt / tt_eft - 1)))
        d["TT max rel vs expected"] = float(np.max(np.abs(tt / tt_exp - 1)))
    cases[name] = d

with warnings.catch_warnings(record=True) as w:
    c = load(); c.set_dark_energy(w=-0.7, dark_energy_model="ppf")
record("set_dark_energy(w=-0.7, ppf)", c, w, do_tt=True)
with warnings.catch_warnings(record=True) as w:
    c = load(); c.set_dark_energy(w=-0.7)
record("set_dark_energy(w=-0.7) [fluid]", c, w)
with warnings.catch_warnings(record=True) as w:
    c = load(); c.set_classes(dark_energy_model="ppf"); c.DarkEnergy.set_params(w=-0.7)
record("set_classes(ppf)+DarkEnergy.set_params(w=-0.7)", c, w)
with warnings.catch_warnings(record=True) as w:
    c = load(); a = np.logspace(-5, 0, 200); c.set_dark_energy_w_a(a, -0.7 * np.ones_like(a), dark_energy_model="ppf")
record("set_dark_energy_w_a(w=-0.7 table, ppf)", c, w)
with warnings.catch_warnings(record=True) as w:
    c = load()
    c = camb.set_params(cp=c, H0=c.H0, ombh2=c.ombh2, omch2=c.omch2, omk=c.omk, mnu=0.0, num_massive_neutrinos=0,
                        nnu=c.N_eff, w=-0.7, dark_energy_model="ppf")
d_sp = {"EFTflag": int(c.EFTCAMB.EFTflag), "H1": H1(c)}
ref = camb.set_params(H0=c.H0, ombh2=c.ombh2, omch2=c.omch2, omk=c.omk, mnu=0.0, num_massive_neutrinos=0,
                      nnu=c.N_eff, w=-0.7, dark_energy_model="ppf")
d_sp["H1 fresh set_params same args"] = H1(ref)
d_sp["rel"] = d_sp["H1"] / d_sp["H1 fresh set_params same args"] - 1
cases["camb.set_params(cp=EFT params, w=-0.7, ppf)"] = d_sp
try:
    with warnings.catch_warnings(record=True) as w:
        c = load()
        p = {"H0": c.H0, "omegabh2": c.ombh2, "omegach2": c.omch2, "tau": 0.06, "w": -0.7,
             "A": 2.1, "ns": 0.965, "mnu": 0.06}
        c = camb.set_params_cosmomc(p, inpars=c, dark_energy_model="ppf")
        c.EFTCAMB.EFTCAMB_feedback_level = 0
    d = {"EFTflag": int(c.EFTCAMB.EFTflag), "H1": H1(c)}
    c.EFTCAMB.EFTflag = 0
    d["H1 same params EFTflag=0"] = H1(c)
    d["rel"] = d["H1"] / d["H1 same params EFTflag=0"] - 1
    cases["set_params_cosmomc(p with w=-0.7, inpars=EFT params)"] = d
except Exception as e:
    cases["set_params_cosmomc(p with w=-0.7, inpars=EFT params)"] = {"error": repr(e)}
out["cases"] = cases
out["TT max rel expected vs EFT"] = float(np.max(np.abs(tt_exp / tt_eft - 1)))
print("@@" + json.dumps(out))
"""


def run(tree: str) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-c", WORKER, TREES + tree, INI, os.path.join(HERE, "..")],
                           capture_output=True, text=True, env=env, preexec_fn=lambda: os.nice(10), timeout=1500)
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@")]
    if not line:
        return {"tree": tree, "error": p.stdout[-2000:] + p.stderr[-3000:]}
    return json.loads(line[0][2:])


if __name__ == "__main__":
    trees = sys.argv[1:] or ["A0", "B1"]
    res = [run(t) for t in trees]
    print(json.dumps(res, indent=1))
    outdir = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
              "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/setde_active_eft")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "result.json"), "w") as f:
        json.dump(res, f, indent=1)

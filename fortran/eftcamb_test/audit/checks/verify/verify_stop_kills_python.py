"""Independent check of finding stop-kills-python via the Python dict API (camb.set_params -> EFTCAMB_initialize in
11_EFTCAMB_python.f90), not camb.read_ini as in the repro.

Each case runs in a child interpreter: set_params + get_results inside try/except Exception; the child prints
CAUGHT <type> if any Python exception arrives and ALIVE after the call. A dead child (no ALIVE) with 'ERROR STOP' in
its output means the Fortran library terminated the interpreter.

Cases include controls where EFTCAMB does use its error mechanism (eft_error -> ValueError, success=.false. ->
CAMBError) to show the contrast in the same code paths.

Usage: python verify_stop_kills_python.py [TREE] [case ...]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REPO = HERE.parents[4]
PAR = REPO / "fortran/eftcamb_test/parameters"
SCRATCH = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
OUT = SCRATCH / "work/verify/stopkill"
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 and "/" in sys.argv[1] else SCRATCH / "trees/A0"
ONLY = [a for a in sys.argv[1:] if "/" not in a]

COSMO = dict(H0=70.0, ombh2=0.0226, omch2=0.112, mnu=0.0, num_massive_neutrinos=0, As=2.1e-9, ns=0.96, tau=0.09,
             WantTensors=False, lmax=2000)
SKIP = {"output_root", "omnuh2"}
if TREE.name == "B1":  # CAMB 1.6.5 set_cosmology fails with numpy 2 on the BBN YHe array: fix YHe explicitly
    COSMO["YHe"] = 0.2454


def eft_keys(ini: str) -> dict:
    """Keys set in the model ini itself (not in base_params.ini)."""
    keys = {}
    for raw in (PAR / ini).read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if "=" in line and not line.upper().startswith("DEFAULT("):
            k, v = (s.strip() for s in line.split("=", 1))
            if k not in SKIP:
                keys[k] = v
    return keys


CHILD = r"""
import json, sys
sys.path.insert(0, sys.argv[1])
cfg = json.loads(sys.argv[2])
import camb
try:
    pars = camb.set_params(**cfg["params"])
    for k, v in cfg.get("attr", {}).items():
        setattr(pars.EFTCAMB, k, v)
    print("SET-OK", flush=True)
    res = camb.get_results(pars)
    print("NO-ERROR TT100=%.6e" % res.get_total_cls(200, CMB_unit="muK")[100, 0], flush=True)
except Exception as e:
    import traceback
    print("TRACE", " | ".join(traceback.format_exc().splitlines()[-6:])[:600], flush=True)
    print("CAUGHT", type(e).__name__, repr(e)[:160], flush=True)
print("ALIVE", flush=True)
"""


def case(model_ini, extra=None, cosmo=None, attr=None, drop=()):
    p = dict(COSMO)
    p.update(cosmo or {})
    p.update({k: v for k, v in eft_keys(model_ini).items() if k not in drop})
    p.update(extra or {})
    p["feedback_level"] = 0
    return {"params": p, "attr": attr or {}}


CASES = {
    # claim (1): Kmouflage + massive neutrinos
    "kmf_nonu_control": case("5_Kmouflage_1.ini"),
    "kmf_mnu": case("5_Kmouflage_1.ini", cosmo=dict(mnu=0.06, num_massive_neutrinos=1)),
    # claim (2): EFT_IC_type outside 1..3, set as a python attribute after initialisation (per-k OMP loop)
    "peft_ic4": case("2_PEFT_Omega_const_1.ini", attr={"EFT_IC_type": 4}, drop=("EFTOmegaExp",)),
    "peft_ic4dict": case("2_PEFT_Omega_const_1.ini", extra={"EFT_IC_type": "4"}, drop=("EFTOmegaExp",)),
    # claim (3)-like: unsupported flag read during initialisation (09_EFTCAMB_main.f90:256)
    "peft_delta_phi": case("2_PEFT_Omega_const_1.ini", extra={"EFTCAMB_evolve_delta_phi": "T"}, drop=("EFTOmegaExp",)),
    # integral() stub of a selectable w_DE parametrization (04p012 exponential_2, flag 11) in designer f(R)
    "desfr_wde11": case("4_Des_F_R_LCDM_1.ini", extra={"EFTwDE": "11"}),
    # control: invalid model flag goes through eft_error -> ValueError, interpreter survives
    "peft_badflag_control": case("2_PEFT_Omega_const_1.ini", extra={"PureEFTmodelOmega": "99"}),
    # control: PEFT nominal through the same python path, and EFT_IC_type=2 (valid) set as attribute
    "peft_control": case("2_PEFT_Omega_const_1.ini", drop=("EFTOmegaExp",)),
    "peft_ic2_control": case("2_PEFT_Omega_const_1.ini", attr={"EFT_IC_type": 2}, drop=("EFTOmegaExp",)),
    # control: scaling cubic nominal
    "scg_control": case("5_scg.ini"),
    # control: scaling cubic omegav < 0.3 -> success=.false. -> CAMBError (same subroutine as the MpiStop)
    "scg_lowomegav_control": case("5_scg.ini", cosmo=dict(omch2=0.36)),
}
# value-dependent: scaling cubic bisection bracket [0, 1.5*4.74e-9] is hard-coded; scan valid parameter values
for A in (-0.3, -0.1, -1.0, -3.0):
    for lam in (154.0, 50.0):
        CASES[f"scg_A{A}_lam{lam}"] = case("5_scg.ini", extra={"Scaling_Cubic_A": str(A), "Scaling_Cubic_lambda": str(lam)})
for b1 in (10.0, 30.0, 300.0, 1000.0):
    CASES[f"scgb_beta1_{b1}"] = case("5_scg.ini", extra={"Scaling_Cubic_beta1": str(b1)})
for b2 in (0.1, 2.0, 5.0):
    CASES[f"scgb_beta2_{b2}"] = case("5_scg.ini", extra={"Scaling_Cubic_beta2": str(b2)})
for oc in (0.30, 0.33):
    CASES[f"scgb_omch2_{oc}"] = case("5_scg.ini", cosmo=dict(omch2=oc))
for oc in (0.05, 0.08, 0.15, 0.25):
    CASES[f"scg_omch2_{oc}"] = case("5_scg.ini", cosmo=dict(omch2=oc))


def run(name, cfg):
    with harness.slot(1):
        r = subprocess.run([sys.executable, "-c", CHILD, str(TREE), json.dumps(cfg)], capture_output=True, text=True,
                           cwd=OUT, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900,
                           preexec_fn=lambda: os.nice(10))
    out = r.stdout + r.stderr
    lines = [ln.strip() for ln in out.splitlines() if ln.strip()]
    res = {"case": name, "rc": r.returncode, "alive": "ALIVE" in r.stdout, "caught": next(
        (ln for ln in lines if ln.startswith("CAUGHT")), None), "noerror": next(
        (ln for ln in lines if ln.startswith("NO-ERROR")), None), "error_stop": "ERROR STOP" in out,
        "tail": lines[-4:]}
    (OUT / f"{TREE.name}_{name}.log").write_text(out)
    return res


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    results = []
    for name, cfg in CASES.items():
        if ONLY and not any(name.startswith(o) for o in ONLY):
            continue
        res = run(name, cfg)
        results.append(res)
        print(f"{name:28s} rc={res['rc']:3d} alive={res['alive']!s:5s} error_stop={res['error_stop']!s:5s} "
              f"caught={res['caught']} noerror={res['noerror']}\n    tail={res['tail']}", flush=True)
    with open(OUT / f"results_{TREE.name}_{'_'.join(ONLY) or 'all'}.json", "w") as f:
        json.dump(results, f, indent=1)


if __name__ == "__main__":
    main()

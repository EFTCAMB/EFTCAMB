"""Independent check for stop-kills-python: use the Python API (set_params + EFTCAMB.initialize_parameters) rather
than read_ini. Each case in its own child process; report whether a Python exception is raised and interpreter lives."""
import json, os, subprocess, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

TREE = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0"
OUT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/stopkill")
PAR = REPO / "fortran/eftcamb_test/parameters"

def eftkeys(ini, **extra):
    d = harness.parse_ini(PAR / ini)
    d = {k: v for k, v in d.items() if k.startswith("EFT") or k in (
        "DesignerEFTmodel", "FullMappingEFTmodel", "PureEFTmodel", "alphaU", "gammaU", "m", "eps2_0", "gammaA")
        or k.startswith("EFTw")}
    d.update({k: str(v) for k, v in extra.items()})
    d["EFTCAMB_write_background"] = "F"
    return d

CASES = [
    ("control: 4_Des_F_R_wCDM_1 as given", eftkeys("4_Des_F_R_wCDM_1.ini"), 0.0),
    ("designer f(R) EFTwDE=7 (Fourier w, valid selector)",
     eftkeys("4_Des_F_R_wCDM_1.ini", EFTwDE=7, EFTw_Fourier_order_cos=1, EFTw_Fourier_order_sin=0,
             EFTw_a0=-0.9, EFTw_a1=0.01), 0.0),
    ("designer f(R) EFTwDE=11 (exponential_2 w)",
     eftkeys("4_Des_F_R_wCDM_1.ini", EFTwDE=11), 0.0),
    ("PEFT Omega const + EFTCAMB_use_background=T via initialize_parameters",
     eftkeys("2_PEFT_Omega_const_1.ini", EFTCAMB_use_background="T"), 0.0),
    ("Kmouflage_1 with mnu=0.06 via set_cosmology", eftkeys("5_Kmouflage_1.ini"), 0.06),
    ("control: unstable PEFT (EFTOmega0=-0.5) -> expected CAMBError",
     eftkeys("2_PEFT_Omega_const_1.ini", EFTOmega0=-0.5), 0.0),
]

CHILD = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
eft = json.loads(sys.argv[2]); mnu = float(sys.argv[3])
try:
    pars = camb.CAMBparams(dark_energy_model="EFTCAMB") if False else camb.set_params(
        H0=67.5, ombh2=0.022, omch2=0.122, mnu=mnu, num_massive_neutrinos=1 if mnu > 0 else 0, As=2e-9, ns=0.965,
        dark_energy_model="EFTCAMB", lmax=1500)
    pars.EFTCAMB.initialize_parameters(pars, eft, print_header=False)
    res = camb.get_results(pars)
    print("NO-ERROR", res.get_cmb_power_spectra(CMB_unit="muK")["total"][200, 0])
except Exception as e:
    print("CAUGHT", type(e).__name__, repr(e)[:160])
print("INTERPRETER-ALIVE")
"""

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for label, eft, mnu in CASES:
        with harness.slot(1):
            r = subprocess.run([sys.executable, "-c", CHILD, TREE, json.dumps(eft), str(mnu)], capture_output=True,
                               text=True, cwd=OUT, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=600)
        tail = [ln.strip() for ln in (r.stdout + r.stderr).splitlines() if ln.strip()][-3:]
        row = dict(case=label, rc=r.returncode, caught="CAUGHT" in r.stdout, alive="INTERPRETER-ALIVE" in r.stdout,
                   tail=tail)
        rows.append(row)
        print(json.dumps(row))
    (OUT / "stopkill_pyapi.json").write_text(json.dumps(rows, indent=1))

if __name__ == "__main__":
    main()

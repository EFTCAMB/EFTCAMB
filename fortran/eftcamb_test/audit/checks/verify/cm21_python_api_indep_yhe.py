"""Independent check for finding 21cm-not-refused-in-python-with-eft.

Different code path from the repro (which used camb.read_ini): build the parameters with the pure Python API
(camb.set_params with EFTCAMB keys + CAMBparams.SourceWindows = [GaussianSourceWindow(source_type="21cm")]), no ini
file involved, and check whether get_results refuses EFTflag != 0 with a 21cm window.
Then run the same model through the Fortran ini driver (camb.run_ini -> CAMB_RunIniFile -> CAMB_RunFromIni) in a
subprocess to confirm that path refuses.

Usage: python cm21_python_api_indep.py [TREE]   (default tree A0)
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

SCR = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit")
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else SCR / "trees" / "A0"
OUT = SCR / "work" / "verify" / "21cm"
OUT.mkdir(parents=True, exist_ok=True)
PARAMS = HERE.parents[2] / "parameters"

COSMO = dict(YHe=0.245, H0=70.0, ombh2=0.0226, omch2=0.112, omk=0.0, tau=0.09, As=2.1e-9, ns=0.96, lmax=600)
EFT_OMEGA = dict(EFTflag=1, PureEFTmodel=1, PureEFTHorndeski=False, EFTwDE=0, PureEFTmodelOmega=1,
                 PureEFTmodelGamma1=0, PureEFTmodelGamma2=0, PureEFTmodelGamma3=0, PureEFTmodelGamma4=0,
                 PureEFTmodelGamma5=0, PureEFTmodelGamma6=0, EFTOmega0=1.0,
                 EFT_ghost_math_stability=False, EFT_mass_math_stability=False, EFT_ghost_stability=True,
                 EFT_gradient_stability=True, EFT_mass_stability=True, EFT_mass_stability_rate=10.0,
                 EFT_additional_priors=True, EFTCAMB_turn_on_time=0.01, EFTCAMB_stability_time=1e-10,
                 EFTCAMB_stability_threshold=0.0)


def python_api_run(eft: dict | None, z: float, sigma_mhz: float = 5.0, kind: str = "21cm") -> dict:
    sys.path.insert(0, str(TREE))
    import camb
    from camb import constants
    from camb.sources import GaussianSourceWindow

    assert Path(camb.__file__).resolve().is_relative_to(TREE.resolve()), camb.__file__
    params = dict(COSMO)
    if eft:
        params.update(eft)
    else:
        params["EFTflag"] = 0
    cp = camb.set_params(**params)
    cp.WantTensors = False
    cp.Want_CMB_lensing = False
    sigma = sigma_mhz / (constants.f_21cm / 1e6) if kind == "21cm" else 0.05
    cp.SourceWindows = [GaussianSourceWindow(redshift=z, source_type=kind, sigma=sigma)]
    cp.SourceTerms.line_basic = True
    cp.SourceTerms.line_distortions = True
    rec = {"EFTflag": int(cp.EFTCAMB.EFTflag), "z": z, "kind": kind}
    try:
        res = camb.get_results(cp)
        cl = res.get_source_cls_dict(raw_cl=True)["W1xW1"]
        rec.update(status="ok", Do21cm=bool(res.Params.Do21cm), cl10=float(cl[10]), cl100=float(cl[100]),
                   cl500=float(cl[500]))
    except Exception as e:  # noqa: BLE001
        rec.update(status=f"exception: {type(e).__name__}: {e}")
    return rec


def ini_driver_run(model: str) -> dict:
    tmp = Path(tempfile.mkdtemp(dir=OUT))
    ini = tmp / f"{model}_21cm.ini"
    ini.write_text(f"DEFAULT({PARAMS / (model + '.ini')})\noutput_root = {tmp}/out\nget_tensor_cls = F\n"
                   "get_transfer = F\nfeedback_level = 1\nEFTCAMB_write_background = F\nDo21cm = T\n"
                   "num_redshiftwindows = 1\nredshift(1) = 1.0\nredshift_kind(1) = 21cm\nredshift_sigma_Mhz(1) = 5\n"
                   "line_basic = T\nline_distortions = T\nl_max_scalar = 600\n")
    code = (f"import sys; sys.path.insert(0, {str(TREE)!r}); import camb; camb.run_ini({str(ini)!r}); "
            "print('RUN_INI_RETURNED_NORMALLY')")
    env = dict(os.environ, OMP_NUM_THREADS="1")
    p = subprocess.run([sys.executable, "-c", code], cwd=tmp, capture_output=True, text=True, timeout=900, env=env)
    txt = p.stdout + p.stderr
    return {"model": model, "rc": p.returncode,
            "refusal_msg": "EFTCAMB is not compatible with 21 cm sources so far" in txt,
            "returned_normally": "RUN_INI_RETURNED_NORMALLY" in txt,
            "tail": txt.strip().splitlines()[-3:]}


if __name__ == "__main__":
    out = {"tree": str(TREE), "python_api": [], "ini_driver": []}
    with harness.slot(1):
        for z in (1.0, 30.0):
            for eft in (None, EFT_OMEGA):
                r = python_api_run(eft, z)
                out["python_api"].append(r)
                print("python API", json.dumps(r))
        for model in ("2_PEFT_Omega_const_1", "1_EFT_GR"):
            r = ini_driver_run(model)
            out["ini_driver"].append(r)
            print("ini driver", json.dumps(r))
    (OUT / ("cm21_python_api_indep_" + TREE.name + "_YHe.json")).write_text(json.dumps(out, indent=1))

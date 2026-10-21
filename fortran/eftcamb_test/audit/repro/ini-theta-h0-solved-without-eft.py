"""Repro: ini keys thetastar / cosmomc_theta (new in CAMB v2) solve H0 before EFTCAMB is initialised.

CAMB_ReadParams -> CAMB_SetH0ForTheta runs with P%EFTCAMB freshly defaulted (EFTflag=0, LCDM background), and the
EFTCAMB model is read only afterwards (CAMB_ReadParamFile / CAMB_RunFromIni).  For EFT models with a non-LCDM
background the H0 written to the params therefore does not reproduce the requested theta.

Expected: theta of the final results == requested theta (brentq rtol 5e-5 on H0, so |rel err| << 1e-4).
Run: /Users/marco/Software/miniforge3/bin/python ini-theta-h0-solved-without-eft.py
"""

import sys
import tempfile
from pathlib import Path

import os  # noqa: E402

os.environ["OMP_NUM_THREADS"] = "2"
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
try:  # hold 2 machine-wide CPU slots of the audit harness for the lifetime of this script, when available
    from harness import slot as _slot

    _SLOT = _slot(2)
    _SLOT.__enter__()
except Exception:  # noqa: BLE001
    pass
TREE = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")
PARAMS = Path(__file__).resolve().parents[2] / "parameters"
sys.path.insert(0, TREE)
import camb  # noqa: E402


def flatten(path: Path) -> dict:
    keys, defaults = {}, []
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.upper().startswith("DEFAULT("):
            defaults.append(flatten(path.parent / line[8:line.rindex(")")].strip()))
        elif "=" in line:
            k, v = line.split("=", 1)
            keys[k.strip()] = v.strip()
    merged = {}
    for d in defaults:
        merged.update(d)
    return merged | keys


THETA = 0.0104110
tmp = Path(tempfile.mkdtemp())
print(f"requested thetastar = {THETA}")
for model in ("1_EFT_GR", "2_PEFT_wDE_wCDM_1", "5_quint_exponential", "5_Kmouflage_1"):
    keys = flatten(PARAMS / f"{model}.ini")
    keys.pop("hubble")
    keys |= {"thetastar": str(THETA), "get_scalar_cls": "F", "get_tensor_cls": "F", "get_transfer": "F",
             "num_redshiftwindows": "0", "feedback_level": "0", "EFTCAMB_write_background": "F"}
    ini = tmp / f"{model}.ini"
    ini.write_text("\n".join(f"{k} = {v}" for k, v in keys.items()) + "\n")
    cp = camb.read_ini(str(ini))
    res = camb.get_background(cp)
    got = res.get_derived_params()["thetastar"] / 100
    # Python path for comparison: same parameters, theta solved through set_cosmology with EFTCAMB initialised
    cp2 = camb.read_ini(str(ini))
    cp2.set_cosmology(H0=None, thetastar=THETA, ombh2=cp2.ombh2, omch2=cp2.omch2, omk=cp2.omk, mnu=None,
                      omnuh2_active=cp2.omnuh2 if cp2.omnuh2 > 0 else 1e-10, nnu=cp2.N_eff, YHe=cp2.YHe,
                      TCMB=cp2.TCMB)
    got2 = camb.get_background(cp2).get_derived_params()["thetastar"] / 100
    print(f"{model:22s} EFTflag={cp.EFTCAMB.EFTflag}  ini: H0={cp.H0:8.4f} thetastar={got:.7f} "
          f"rel.err={got / THETA - 1:+.2e} | set_cosmology: H0={cp2.H0:8.4f} rel.err={got2 / THETA - 1:+.2e}")

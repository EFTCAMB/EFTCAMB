"""Extra controls for finding eft-par-cache-realloc-crash.

1) B0 (plain upstream CAMB v2): results.get_BAO(z) with default params works (expected API behaviour).
2) A0: other API paths that feed results.Params back into Fortran SetParams
   (calc_power_spectra(results.Params), get_results(results.Params) with full Cls), and whether
   the abort can be caught from Python (try/except around get_BAO).
Each case in a subprocess. Usage: python verify_par_cache_realloc_extra.py
"""

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
OUT = ROOT / "work" / "verify" / "par_cache_realloc" / "extra"

CHILD = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
mode, ini = sys.argv[2], sys.argv[3]
z = [0.5, 1.0]
if mode == "B0_get_BAO_default":
    cp = camb.set_params(H0=70, ombh2=0.0226, omch2=0.112, mnu=0.06)
    res = camb.get_background(cp)
    print("RESULT", json.dumps([res.hubble_parameter(z).tolist(), res.get_BAO(z)[:, 1].tolist()]), flush=True)
    sys.exit(0)
cp = camb.read_ini(ini, no_validate=True)
if mode.endswith("_GR"):
    cp.EFTCAMB.EFTflag = 0
    mode = mode[:-3]
print("EFTFLAG", int(cp.EFTCAMB.EFTflag), flush=True)
if mode == "calc_power_spectra_resParams":
    res = camb.get_results(cp)
    ref = float(res.get_cmb_power_spectra(CMB_unit="muK")["total"][100, 0])
    res.calc_power_spectra(res.Params)
    print("RESULT", json.dumps([ref, float(res.get_cmb_power_spectra(CMB_unit="muK")["total"][100, 0])]), flush=True)
elif mode == "try_except_get_BAO":
    res = camb.get_background(cp)
    try:
        out = res.get_BAO(z)[:, 1].tolist()
    except Exception as e:  # noqa: BLE001
        out = f"caught {type(e).__name__}: {e}"
    print("RESULT", json.dumps(out), flush=True)
"""


def run(tree, mode, ini):
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", CHILD, str(tree), mode, str(ini)], capture_output=True, text=True,
                       env=env, cwd=str(OUT))
    res = [ln[7:] for ln in p.stdout.splitlines() if ln.startswith("RESULT")]
    err = [ln for ln in (p.stdout + p.stderr).splitlines() if "Fortran runtime error" in ln or "At line" in ln
           or "Error" in ln]
    return {"tree": Path(tree).name, "mode": mode, "rc": p.returncode, "result": json.loads(res[0]) if res else None,
            "err": err[-2:]}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ini = OUT / "fR.ini"
    ini.write_text(f"DEFAULT({PARAMS / '4_Des_F_R_LCDM_1.ini'})\noutput_root = {OUT}/fR\n"
                   "EFTCAMB_write_background = F\nfeedback_level = 0\n")
    cases = [(ROOT / "trees" / "B0", "B0_get_BAO_default"),
             (ROOT / "trees" / "A0", "try_except_get_BAO"),
             (ROOT / "trees" / "A0", "calc_power_spectra_resParams"),
             (ROOT / "work" / "verify" / "fixmin_realloc", "calc_power_spectra_resParams"),
             (ROOT / "trees" / "A0", "calc_power_spectra_resParams_GR"),
             (ROOT / "trees" / "B1", "calc_power_spectra_resParams_GR")]
    allres = []
    for tree, mode in cases:
        r = run(tree, mode, ini)
        allres.append(r)
        print(json.dumps(r).replace(str(ROOT), "$S"), flush=True)
    (OUT / "extra.json").write_text(json.dumps(allres, indent=1))


if __name__ == "__main__":
    main()

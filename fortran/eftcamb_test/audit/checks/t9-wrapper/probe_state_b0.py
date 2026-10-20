"""Track 9 (v) control: dependence of a run on a previous run (different lmax/settings) in the same process.

Reference R = read_ini(DEFAULT(<model>.ini) + l_max_scalar = 1200, get_tensor_cls = F), computed fresh and after a
GR run from a fresh CAMBparams at lmax 1000 (lens_potential_accuracy 0).  Cases:
  B0/GR : upstream CAMB, model 1_EFT_GR.ini (EFT keys ignored by B0)
  A0/GR : code under audit, same GR ini
  A0/EFT: code under audit, 2_PEFT_Omega_power_1.ini (the probe_errors reference)
"""

import json
import sys

import numpy as np

from t9lib import PARAMS, WORK, run_snippet

CODE = r"""
ini = WORK / 'state_{tag}.ini'
ini.write_text('DEFAULT({model})\nl_max_scalar = 1200\nget_tensor_cls = F\nEFTCAMB_write_background = F\n')
if {after}:
    p = camb.CAMBparams()
    p.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12)
    p.set_for_lmax(1000, lens_potential_accuracy=0)
    camb.get_results(p)
r = camb.get_results(camb.read_ini(str(ini)))
np.save(str(WORK / 'state_{tag}_{after}.npy'), r.get_cmb_power_spectra(raw_cl=True, CMB_unit='muK')['total'])
emit(done=True)
"""

CASES = {"B0_GR": ("B0", "1_EFT_GR.ini"), "A0_GR": ("A0", "1_EFT_GR.ini"), "A0_EFT": ("A0", "2_PEFT_Omega_power_1.ini")}

if __name__ == "__main__":
    out = {}
    for tag in sys.argv[1:] or list(CASES):
        tree, model = CASES[tag]
        for after in (False, True):
            r = run_snippet(CODE.format(tag=tag, model=PARAMS / model, after=after), tree_name=tree, threads=2,
                            timeout=600)
            assert r["status"] == "ok", r["tail"]
        a = np.load(WORK / f"state_{tag}_False.npy")
        b = np.load(WORK / f"state_{tag}_True.npy")
        m = np.abs(a) > 0
        out[tag] = {"identical": bool(np.array_equal(a, b)), "maxrel": float(np.max(np.abs(a[m] - b[m]) / np.abs(a[m])))}
        print(tag, out[tag], flush=True)

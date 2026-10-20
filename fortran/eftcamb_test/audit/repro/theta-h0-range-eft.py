"""Repro (track t8c-features): H0 search range for theta with EFTCAMB.

camb/model.py set_cosmology (EFTCAMB MOD) raises the lower end of theta_H0_range to 1.1*100*sqrt(Omega_m h^2)
(Omega_m(H0_min) = 0.826) for EVERY model, including EFTflag = 0 and omk != 0, while CAMBparams.set_H0_for_theta (the
public method it wraps, and the Fortran CAMB_SetH0ForTheta used by the ini keys) keeps (10, 100).

 1. Designer f(R) (4_Des_F_R_LCDM_1): cp.set_H0_for_theta(thetastar) with the default range evaluates H0 = 10, where
    Omega_DE < 0 and the EFT model is unstable; the CAMBError "EFTCAMB: theory unstable" aborts the root search
    although a stable solution exists in the range (found by set_cosmology(thetastar=...) and by an explicit range).
 2. GR (EFTflag = 0), Omega_m(H0=58) = 0.96: set_cosmology(thetastar=...) reports "No solution" in A0 while plain CAMB
    v2 (B0) recovers H0 = 58.
Expected: same behaviour as CAMB for GR, and set_H0_for_theta usable with EFT models.
Run: /Users/marco/Software/miniforge3/bin/python theta-h0-range-eft.py
"""

import json
import subprocess
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
ROOT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
            "scratchpad/audit/trees")
PARAMS = Path(__file__).resolve().parents[2] / "parameters"

CASE1 = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
ini = sys.argv[2]
out = {}
theta = 0.010411
for tag, kw in (("set_H0_for_theta default range", {}), ("set_H0_for_theta range (55,100)", {"theta_H0_range": (55, 100)})):
    cp = camb.read_ini(ini)
    try:
        cp.set_H0_for_theta(theta, **kw)
        r = camb.get_background(cp)
        out[tag] = "H0=%.4f thetastar=%.7f" % (cp.H0, r.get_derived_params()["thetastar"] / 100)
    except Exception as e:
        out[tag] = "%s: %s" % (type(e).__name__, str(e).splitlines()[-1])
cp = camb.read_ini(ini)
cp.set_cosmology(H0=None, thetastar=theta, ombh2=cp.ombh2, omch2=cp.omch2, omk=cp.omk, mnu=None,
                 omnuh2_active=cp.omnuh2, num_massive_neutrinos=max(int(cp.num_nu_massive), 1), nnu=cp.N_eff,
                 YHe=cp.YHe, TCMB=cp.TCMB)
r = camb.get_background(cp)
out["set_cosmology(thetastar)"] = "H0=%.4f thetastar=%.7f" % (cp.H0, r.get_derived_params()["thetastar"] / 100)
print(json.dumps(out))
"""

CASE2 = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
p = camb.CAMBparams()
p.set_cosmology(H0=58, ombh2=0.0226, omch2=0.3, mnu=0.06)
om = (p.ombh2 + p.omch2 + p.omnuh2) / (p.H0 / 100) ** 2
target = camb.get_background(p).get_derived_params()["thetastar"] / 100
q = camb.CAMBparams()
try:
    q.set_cosmology(thetastar=target, H0=None, ombh2=0.0226, omch2=0.3, mnu=0.06)
    res = "H0=%.4f" % q.H0
except Exception as e:
    res = "%s: %s" % (type(e).__name__, e)
print(json.dumps({"Omega_m(H0=58)": round(om, 3), "target thetastar": target, "set_cosmology(thetastar)": res}))
"""

tmp = Path(tempfile.mkdtemp())
ini = tmp / "fr.ini"
ini.write_text(f"DEFAULT({PARAMS / '4_Des_F_R_LCDM_1.ini'})\nget_scalar_cls = F\nget_tensor_cls = F\nget_transfer = F\n"
               "num_redshiftwindows = 0\nfeedback_level = 0\nEFTCAMB_write_background = F\n")


def run(code, *args):
    p = subprocess.run([sys.executable, "-c", code, *map(str, args)], capture_output=True, text=True, cwd=tmp,
                       env={"OMP_NUM_THREADS": "2", "PATH": "/usr/bin:/bin"}, preexec_fn=lambda: os.nice(10))
    lines = [x for x in p.stdout.splitlines() if x.startswith("{")]
    return json.loads(lines[-1]) if lines else {"rc": p.returncode, "stderr": p.stderr[-400:]}


print("Case 1: designer f(R) LCDM (A0), requested thetastar = 0.010411")
for k, v in run(CASE1, ROOT / "A0", ini).items():
    print(f"   {k:34s}: {v}")
print("Case 2: GR, high Omega_m")
for t in ("A0", "B0"):
    print(f"   {t}: {run(CASE2, ROOT / t)}")

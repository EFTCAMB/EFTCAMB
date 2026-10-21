"""Verifier: compare the model's own analytic second derivatives (EFT_mass_stability=F, so
compute_additional_derivs is not called) with what EFTTestStability leaves in the cache when
EFT_mass_stability=T (after compute_additional_derivs). Also FD of analytic first derivatives.
Usage: python addder_analytic_vs_overwrite.py TREE"""
import os, sys, json, subprocess, tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness
TREE = sys.argv[1]
PAR = Path(__file__).resolve().parents[3] / "parameters"
CHILD = r'''
import sys, json, ctypes
sys.path.insert(0, sys.argv[1])
import camb
from camb import eftcamb
from camb.baseconfig import camblib
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.EFTCAMB.EFT_mass_stability = False
res = camb.get_results(pars); P = res.Params
f = camblib.__eftcamb_stability_MOD_eftteststability; f.restype = ctypes.c_int
def cache(a, mass):
    P.EFTCAMB.EFT_mass_stability = mass
    c = eftcamb.EFTCAMB_timestep_cache()
    ok = f(ctypes.byref(ctypes.c_double(a)), ctypes.byref(ctypes.c_double(10.0)), ctypes.byref(P.EFTCAMB),
           ctypes.byref(P.EFTCAMB_parameter_cache), ctypes.byref(c))
    return c, bool(ok)
out = []
for a in [0.0101, 0.1, 0.5]:
    cF, _ = cache(a, False); cT, okT = cache(a, True)
    eps = 1e-5; h = 2*a*eps
    cp, _ = cache(a*(1+eps), False); cm, _ = cache(a*(1-eps), False)
    r = dict(a=a, stable_mass=okT, mu1_over_H2=cT.EFT_mu1/cT.adotoa**2)
    for k, kp in (("EFTGamma1PP", "EFTGamma1P"), ("EFTGamma2PP", "EFTGamma2P"), ("EFTGamma3PP", "EFTGamma3P")):
        r[k] = [getattr(cF, k), getattr(cT, k), (getattr(cp, kp) - getattr(cm, kp))/h]
    fd = a*cF.adotoa*(cp.EFTcdot - cm.EFTcdot)/h - 2*cF.adotoa*cF.EFTcdot
    r["EFTcdotdot"] = [cF.EFTcdotdot, cT.EFTcdotdot, fd, -2*cF.adotoa*cF.EFTcdot]
    out.append(r)
print("JSON" + json.dumps(out))
'''
tmp = Path(tempfile.mkdtemp(prefix="addderV_", dir=os.environ.get("OUTDIR")))
for ini in ["3_RPH_KT_power_1.ini", "3_RPH_MK_power_1.ini", "2_PEFT_Omega_const_1.ini"]:
    w = tmp / ini
    w.write_text(f"DEFAULT({PAR/ini})\nfeedback_level = 0\nEFTCAMB_write_background = F\nget_tensor_cls = F\n"
                 f"do_lensing = F\noutput_root = {tmp/ini}\n")
    with harness.slot(1):
        r = subprocess.run([sys.executable, "-c", CHILD, TREE, str(w)], capture_output=True, text=True, cwd=tmp,
                           env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900)
    js = [l for l in r.stdout.splitlines() if l.startswith("JSON")]
    print(f"== {ini} rc={r.returncode}")
    if not js:
        print(r.stdout[-1000:], r.stderr[-1500:]); continue
    for row in json.loads(js[-1][4:]):
        print(f" a={row['a']}  mass-stable={row['stable_mass']}  mu1/H^2(used)={row['mu1_over_H2']:+.4e}")
        for k in ("EFTGamma1PP", "EFTGamma2PP", "EFTGamma3PP", "EFTcdotdot"):
            v = row[k]
            print(f"   {k:12s} analytic={v[0]:+.6e}  after_additional_derivs={v[1]:+.6e}  FD(analytic 1st)={v[2]:+.6e}"
                  + (f"  -2Hcdot={v[3]:+.6e}" if len(v) > 3 else ""))

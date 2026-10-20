"""Independent check for theta-solve-aborts-unstable-bracket-end: use set_cosmology's own theta_H0_range
argument (no hand-written brentq) and the cosmomc_theta from CAMBdata after a full get_results; scan H0 for
stability boundary; run SCG direct H0 in subprocess (A0 and B1)."""
import json, os, subprocess, sys
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
TREES = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/"
PARAMS = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"
W = r"""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1]); import camb
ini, mode = sys.argv[2], sys.argv[3]
cp = camb.read_ini(ini, no_validate=True); camb.set_feedback_level(0)
try: cp.EFTCAMB.EFTCAMB_feedback_level = 0
except Exception: pass
r = {}
if mode == "direct":
    cp.H0 = float(sys.argv[4]); print("start", flush=True)
    try: camb.get_background(cp, no_thermo=True); r["res"] = "ok"
    except camb.CAMBError as e: r["res"] = str(e).splitlines()[-1]
else:
    th = camb.get_background(cp).cosmomc_theta(); r["theta"] = th
    kw = dict(ombh2=cp.ombh2, omch2=cp.omch2, omk=cp.omk, mnu=None, nnu=cp.N_eff, TCMB=cp.TCMB, YHe=cp.YHe,
              omnuh2_active=cp.omnuh2, num_massive_neutrinos=max(1, int(cp.num_nu_massive)) if cp.omnuh2 > 0 else 0)
    scan = {}
    for H in [72, 75, 80, 85, 90, 95, 100]:
        cp.H0 = H
        try: camb.get_background(cp, no_thermo=True); scan[H] = "ok"
        except camb.CAMBError as e: scan[H] = str(e).splitlines()[-1]
    r["scan"] = scan
    for hi in (100, 74):
        try:
            cp.set_cosmology(cosmomc_theta=th, theta_H0_range=(10, hi), **kw); r["hi%d" % hi] = cp.H0
        except Exception as e: r["hi%d" % hi] = repr(e).splitlines()[-1][:90]
print("@@" + json.dumps(r))
"""
env = dict(os.environ, OMP_NUM_THREADS="1"); env.pop("PYTHONPATH", None)
def run(args):
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-c", W] + args, capture_output=True, text=True, env=env,
                           preexec_fn=lambda: os.nice(10), timeout=600)
    l = [x for x in p.stdout.splitlines() if x.startswith("@@")]
    return json.loads(l[0][2:]) if l else {"died_rc": p.returncode, "tail": (p.stdout + p.stderr).strip().splitlines()[-3:]}
for m in ["2_PEFT_Omega_power_2", "3_RPH_KT_OmegaDE_1", "5_quint_inverse_power"]:
    print(m, "A0", run([TREES + "A0", f"{PARAMS}/{m}.ini", "theta"]))
for T in ["A0", "B1"]:
    for H in ["70", "95", "100"]:
        print("5_scg", T, H, run([TREES + T, f"{PARAMS}/5_scg.ini", "direct", H]))

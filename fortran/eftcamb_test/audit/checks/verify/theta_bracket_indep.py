"""Independent verification of finding 'theta-solve-aborts-unstable-bracket-end'.

Method (different from the repro):
  1. 'scan': for each model, load the ini with read_ini (no set_cosmology rebuild), and scan H0 on a grid between
     1.1*100*sqrt(omega_m) and 100, calling camb.get_background(no_thermo=True) on a deep copy for every point;
     record ok / CAMBError message / process death (MpiStop).  This maps the valid H0 region and checks that
     theta(H0) is monotonic there, so that the target theta has a unique root inside the valid region.
  2. 'api': public API only: set_cosmology(cosmomc_theta=theta_ini) with the default theta_H0_range, the same
     with thetastar (cosmomc_approx=False path), and set_cosmology(cosmomc_theta, theta_H0_range=<valid sub-range>)
     (the documented workaround), compared with H0_ini.
  3. controls: B0 (upstream CAMB, LCDM same densities) and B1 (pre-rebase EFTCAMB) for the same calls.

Usage: python theta_bracket_indep.py <outdir> [models...]
"""

import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import harness  # noqa: E402

TREES = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees")
PARAMS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "parameters"))
CHECKS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

WORKER = r"""
import sys, json, copy
sys.path.insert(0, sys.argv[1])
sys.path.insert(0, sys.argv[2])
import numpy as np
import camb
import harness
task, ini = sys.argv[3], sys.argv[4]
H0s = json.loads(sys.argv[5]) if len(sys.argv) > 5 else []
camb.set_feedback_level(0)
keys = harness.parse_ini(ini)
def fresh():
    # a fresh parameter object read from the ini for every evaluation (no state carried between trials)
    p = camb.read_ini(ini, no_validate=True)
    if hasattr(p, "EFTCAMB"):
        if p.EFTCAMB.EFTflag == 0 and int(keys.get("EFTflag", "0")) != 0:
            p.EFTCAMB.initialize_parameters(p, keys, print_header=False)   # B1 fallback
        p.EFTCAMB.EFTCAMB_feedback_level = 0
    return p
cp = fresh()
has_eft = hasattr(cp, "EFTCAMB")
def theta_at(H0):
    p = fresh()
    p.H0 = H0
    return camb.get_background(p, no_thermo=True).cosmomc_theta()
if task == "scan":
    for H0 in H0s:
        print("@@TRY %r" % H0, flush=True)
        try:
            th = theta_at(H0)
            print("@@RES " + json.dumps([H0, "ok", th]), flush=True)
        except camb.CAMBError as e:
            print("@@RES " + json.dumps([H0, "err", str(e).splitlines()[-1]]), flush=True)
elif task == "api":
    H0_ini = cp.H0
    out = {"H0_ini": H0_ini, "EFTflag": int(cp.EFTCAMB.EFTflag) if has_eft else None}
    theta = theta_at(H0_ini)
    out["theta_ini"] = theta
    # theta* (cosmomc_approx=False) target
    bg = camb.get_background(fresh())
    out["thetastar_ini"] = bg.get_derived_params()["thetastar"] / 100.0
    kw = dict(ombh2=cp.ombh2, omch2=cp.omch2, omk=cp.omk, mnu=None, nnu=cp.N_eff, TCMB=cp.TCMB, YHe=cp.YHe,
              omnuh2_active=cp.omnuh2, num_massive_neutrinos=int(cp.num_nu_massive) if cp.omnuh2 > 0 else 0)
    import inspect
    if "omnuh2_active" not in inspect.signature(cp.set_cosmology).parameters:   # CAMB 1.6.5 (B1): mnu only
        from camb import constants
        kw.pop("omnuh2_active")
        nfac = constants.neutrino_mass_fac * (constants.COBE_CMBTemp / cp.TCMB) ** 3
        kw["mnu"] = cp.omnuh2 * nfac / (constants.default_nnu / 3) ** 0.75
    # consistency of rebuilt parameters: set_cosmology(H0=H0_ini) gives the same theta as the ini
    p = fresh()
    p.set_cosmology(H0=H0_ini, **kw)
    out["theta_rebuild_rel"] = camb.get_background(p, no_thermo=True).cosmomc_theta() / theta - 1
    rng = json.loads(sys.argv[5]) if len(sys.argv) > 5 else None
    cases = [("cosmomc_theta_default", dict(cosmomc_theta=theta)),
             ("thetastar_default", dict(thetastar=out["thetastar_ini"]))]
    if rng:
        cases.append(("cosmomc_theta_range%s" % rng, dict(cosmomc_theta=theta, theta_H0_range=tuple(rng))))
    for name, extra in cases:
        print("@@TRY " + name, flush=True)
        p = fresh()
        try:
            p.set_cosmology(**kw, **extra)
            out[name] = p.H0
        except Exception as e:
            out[name] = type(e).__name__ + ": " + str(e).splitlines()[-1]
    print("@@API " + json.dumps(out), flush=True)
"""


def _env():
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    return env


def _call(args):
    with harness.slot(1):
        return subprocess.run([sys.executable, "-c", WORKER, *args], capture_output=True, text=True, env=_env(),
                              preexec_fn=lambda: os.nice(10), timeout=600)


def scan(tree, ini, H0s):
    """Scan theta(H0); restart after a process death (MpiStop) with the remaining points."""
    res = []
    todo = list(H0s)
    while todo:
        p = _call([tree, CHECKS, "scan", ini, json.dumps(todo)])
        done = [json.loads(ln[6:]) for ln in p.stdout.splitlines() if ln.startswith("@@RES ")]
        res += done
        todo = todo[len(done):]
        if todo:  # died on todo[0]
            tail = [ln.strip() for ln in (p.stdout + p.stderr).splitlines() if ln.strip() and not ln.startswith("@@")]
            res.append([todo[0], "DIED rc=%d" % p.returncode, " | ".join(tail[-3:])])
            todo = todo[1:]
    return res


def api(tree, ini, rng=None):
    args = [tree, CHECKS, "api", ini] + ([json.dumps(rng)] if rng else [])
    p = _call(args)
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@API ")]
    if line:
        return json.loads(line[0][6:])
    tries = [ln for ln in p.stdout.splitlines() if ln.startswith("@@TRY")]
    tail = [ln.strip() for ln in (p.stdout + p.stderr).splitlines() if ln.strip() and not ln.startswith("@@")]
    return {"DIED": "rc=%d during %s: %s" % (p.returncode, tries[-1:] if tries else "?", " | ".join(tail[-3:]))}


if __name__ == "__main__":
    outdir = sys.argv[1]
    mode = sys.argv[2]
    models = sys.argv[3:]
    os.makedirs(outdir, exist_ok=True)
    if mode == "scan":
        import numpy as np
        out = {}
        def job(m):
            ini = os.path.join(PARAMS, m + ".ini")
            keys = harness.parse_ini(ini)
            om = float(keys["ombh2"]) + float(keys["omch2"]) + float(keys.get("omnuh2", 0))
            lo = 1.1 * 100 * np.sqrt(om)
            grid = [round(float(x), 4) for x in np.concatenate([[lo], np.arange(45, 100, 2.5), [100.0]])]
            return m, scan(os.path.join(TREES, "A0"), ini, grid)
        with ThreadPoolExecutor(4) as ex:
            for m, r in ex.map(job, models):
                out[m] = r
                print(m)
                for row in r:
                    print("   ", row)
        json.dump(out, open(os.path.join(outdir, "scan_A0.json"), "w"), indent=1)
    elif mode == "api":
        tree = sys.argv[3]
        rng = json.loads(sys.argv[4])
        models = sys.argv[5:]
        out = {}
        def job(m):
            ini = m if m.endswith(".ini") else os.path.join(PARAMS, m + ".ini")
            return m, api(os.path.join(TREES, tree), ini, rng.get(m))
        with ThreadPoolExecutor(4) as ex:
            for m, r in ex.map(job, models):
                out[m] = r
                print(tree, m, json.dumps(r, indent=None))
        json.dump(out, open(os.path.join(outdir, "api_%s.json" % tree), "w"), indent=1)

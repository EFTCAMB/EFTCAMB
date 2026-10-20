"""Independent check of finding eft-pickle-repr-loses-model with EFT-sensitive observables and several models.

For each model ini (read with camb.read_ini on the chosen tree), in separate subprocesses:
  ref      : get_results on the original params -> sigma8, H(z=0.5,1), TT(ell=2,10) (muK^2)
  fields   : pickle round trip; compare all python-mapped EFTCAMB fields (must be equal) and model_name()
  crash    : pickle round trip -> CAMBdata.calc_background_no_thermo (no thermo, no get_results)
  reinit   : pickle round trip, then EFTCAMB.initialize_parameters(flattened ini) -> get_results; must equal ref
  forkpool : ProcessPoolExecutor with the *fork* start method (arguments are still pickled) -> sigma8
GR control: EFTflag forced to 0 after read_ini, pickle -> get_results must equal the unpickled GR run.

Usage: python verify_pickle_model_obs.py [--tree PATH] [--models a,b,c]
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
REPO = Path(__file__).resolve().parents[5]
PARAMS = REPO / "fortran/eftcamb_test/parameters"
WD = Path(f"{ROOT}/work/verify/pickle_obs")
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))

CHILD = r"""
import sys, json, pickle
tree, mode, ini, flat_json = sys.argv[1:5]
sys.path.insert(0, tree)
import camb
flat = json.loads(flat_json)

def obs(p):
    r = camb.get_results(p)
    tt = r.get_cmb_power_spectra(p, CMB_unit="muK", raw_cl=True)["total"][:, 0]
    return {"sigma8": float(r.get_sigma8_0()), "H": r.hubble_parameter([0.5, 1.0]).tolist(),
            "TT2": float(tt[2]), "TT10": float(tt[10])}

def work(p):
    return obs(p)

if __name__ == "__main__":
    p = camb.read_ini(ini, no_validate=True)
    p.WantTransfer = True
    p.set_matter_power(redshifts=[0.0], kmax=2.0)
    out = {}
    if mode == "ref":
        out = obs(p)
    elif mode == "gr":
        p.EFTCAMB.EFTflag = 0
        a = obs(p)
        q = pickle.loads(pickle.dumps(p))
        out = {"gr_ref": a, "gr_pickle": obs(q)}
    elif mode == "fields":
        q = pickle.loads(pickle.dumps(p))
        names = [n for n, _ in p.EFTCAMB.get_all_fields()]
        diff = [n for n in names if getattr(p.EFTCAMB, n) != getattr(q.EFTCAMB, n)]
        out = {"n_fields": len(names), "differing_fields": diff, "model_name_orig": p.EFTCAMB.model_name(),
               "model_name_restored": q.EFTCAMB.model_name(), "EFTflag_restored": q.EFTCAMB.EFTflag}
    elif mode == "crash":
        q = pickle.loads(pickle.dumps(p))
        print("BEFORE calc_background_no_thermo", flush=True)
        d = camb.CAMBdata()
        d.calc_background_no_thermo(q)
        out = {"H": d.hubble_parameter([0.5, 1.0]).tolist()}
    elif mode == "reinit":
        q = pickle.loads(pickle.dumps(p))
        q.EFTCAMB.initialize_parameters(q, flat, print_header=False)
        out = obs(q)
    elif mode == "forkpool":
        import concurrent.futures as cf, multiprocessing as mp
        try:
            with cf.ProcessPoolExecutor(1, mp_context=mp.get_context("fork")) as ex:
                out = ex.submit(work, p).result(timeout=300)
        except Exception as e:
            out = {"pool_exception": type(e).__name__}
    print("JSON " + json.dumps(out), flush=True)
"""


def run(tree, mode, ini, flat):
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    child = WD / "child_obs.py"
    p = subprocess.run([sys.executable, str(child), tree, mode, str(ini), json.dumps(flat)], capture_output=True,
                       text=True, env=env, cwd=WD, timeout=900)
    res = None
    for ln in p.stdout.splitlines():
        if ln.startswith("JSON "):
            res = json.loads(ln[5:])
    return p.returncode, res, p.stderr.strip()[-200:]


def rel(a, b):
    return abs(a - b) / max(abs(a), 1e-300)


def main():
    import harness

    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=f"{ROOT}/trees/A0")
    ap.add_argument("--models", default="2_PEFT_wDE_wCDM_1,4_Des_mc5e_CPL_1,5_hdsk_TG_1")
    a = ap.parse_args()
    WD.mkdir(parents=True, exist_ok=True)
    (WD / "child_obs.py").write_text(CHILD)
    summary = {}
    for m in a.models.split(","):
        ini = WD / f"{m}.ini"
        ini.write_text(f"DEFAULT({PARAMS / (m + '.ini')})\noutput_root = {WD}/{m}\n"
                       "EFTCAMB_write_background = F\nfeedback_level = 0\nEFTCAMB_feedback_level = 0\n")
        flat = harness.parse_ini(ini)
        flat = {k: v for k, v in flat.items() if k != "output_root"}
        res = {}
        for mode in ("ref", "fields", "crash", "reinit", "forkpool", "gr"):
            rc, r, err = run(a.tree, mode, ini, flat)
            res[mode] = {"rc": rc, "res": r, "err": err if rc else ""}
            st = "OK" if rc == 0 else (f"signal {-rc}" if rc < 0 else f"rc={rc}")
            print(f"{m:22s} {mode:9s} {st:10s} {json.dumps(r)} {res[mode]['err']}", flush=True)
        ref, ri = res["ref"]["res"], res["reinit"]["res"]
        if ref and ri:
            print(f"{m:22s} reinit-vs-ref max rel diff: "
                  f"{max(rel(ref[k], ri[k]) for k in ('sigma8', 'TT2', 'TT10')):.3e}")
        g = res["gr"]["res"]
        if ref and g:
            print(f"{m:22s} EFT effect (ref vs GR): sigma8 {rel(g['gr_ref']['sigma8'], ref['sigma8']):.3e}, "
                  f"TT2 {rel(g['gr_ref']['TT2'], ref['TT2']):.3e}; GR pickle vs GR: "
                  f"{max(rel(g['gr_ref'][k], g['gr_pickle'][k]) for k in ('sigma8', 'TT2', 'TT10')):.3e}")
        summary[m] = res
    (WD / f"summary_{Path(a.tree).name}.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

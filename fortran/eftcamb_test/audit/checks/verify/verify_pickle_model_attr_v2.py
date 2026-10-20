"""Independent check (reproduce + attribute) of finding eft-pickle-repr-loses-model.

Different path from the repro: parameters are built with camb.set_params (or read_ini with EFTflag reset to 0)
plus EFTCAMB.initialize_parameters from python; results via camb.get_results (full SetParams) and
CAMBdata.calc_background_no_thermo; the pickle round trip that multiprocessing does (ProcessPoolExecutor, spawn).
Also inspects the raw memory of TEFTCAMB beyond the ctypes mapping (A0 layout: 328 bytes, ctypes maps 168),
GR controls and a re-initialisation workaround.

Usage: python verify_pickle_model_attr_v2.py [--tree PATH] [--ini INI] [--modes m1,m2]
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
REPO = Path(__file__).resolve().parents[5]
FR_INI = REPO / "fortran/eftcamb_test/parameters/4_Des_F_R_LCDM_1.ini"
WD = Path(f"{ROOT}/work/verify/pickle_model_attr")
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))

CHILD = r"""
import sys, json, ctypes, pickle, copy
tree, mode, eft_json = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, tree)
import camb
eft = json.loads(eft_json)
out = {}

def make(eftflag=True):
    if len(sys.argv) > 4:  # build from read_ini (B1: set_params broken with this numpy); EFT from python
        p = camb.read_ini(sys.argv[4], no_validate=True)
        p.EFTCAMB.EFTflag = 0
    else:
        p = camb.set_params(H0=70, ombh2=0.0226, omch2=0.112, mnu=0.06, As=2.1e-9, ns=0.96)
    if eftflag:
        p.EFTCAMB.initialize_parameters(p, eft, print_header=False)
    return p

def raw_tail(e):
    n = ctypes.sizeof(type(e))
    buf = (ctypes.c_char * 328).from_address(ctypes.addressof(e)).raw
    words = [int.from_bytes(buf[i:i + 8], "little") for i in range(168, 328, 8)]
    return n, [i for i, w in enumerate(words) if w != 0]

def hub(r):
    return r.hubble_parameter([0.5, 1.0]).tolist()

if mode == "state":
    p = make()
    st = p.__getstate__()["EFTCAMB"]
    n, nz = raw_tail(p.EFTCAMB)
    q = pickle.loads(pickle.dumps(p))
    n2, nz2 = raw_tail(q.EFTCAMB)
    d = copy.deepcopy(p)
    n3, nz3 = raw_tail(d.EFTCAMB)
    out.update(ctypes_sizeof=n, eft_state_type=type(st).__name__, n_state_fields=len(st.__getstate__()),
               orig_nonzero_words_after_168=nz, pickled_nonzero_words_after_168=nz2,
               deepcopy_nonzero_words_after_168=nz3, pickled_EFTflag=q.EFTCAMB.EFTflag)
elif mode == "ref":
    r = camb.get_results(make())
    out.update(H=hub(r), age=r.get_derived_params()["age"])
elif mode == "gr_pickle":
    p = make(eftflag=False)
    out.update(H_ref=hub(camb.get_background(p)), H_pickle=hub(camb.get_background(pickle.loads(pickle.dumps(p)))),
               H_deepcopy=hub(camb.get_background(copy.deepcopy(p))))
elif mode == "deepcopy_results":
    q = copy.deepcopy(make())
    print("BEFORE get_results EFTflag=%d" % q.EFTCAMB.EFTflag, flush=True)
    r = camb.get_results(q)
    out.update(H=hub(r))
elif mode == "set_params_cambdata":
    q = pickle.loads(pickle.dumps(make()))
    print("BEFORE CAMBdata.calc_background_no_thermo", flush=True)
    d = camb.CAMBdata()
    d.calc_background_no_thermo(q)
    out.update(H=hub(d))
elif mode == "workaround":
    q = pickle.loads(pickle.dumps(make()))
    q.EFTCAMB.initialize_parameters(q, eft, print_header=False)
    r = camb.get_results(q)
    out.update(H=hub(r), age=r.get_derived_params()["age"])
elif mode == "pool":
    import concurrent.futures as cf, multiprocessing as mp
    import _verify_pool_worker_attr as w
    p = make()
    try:
        with cf.ProcessPoolExecutor(1, mp_context=mp.get_context("spawn"), initializer=w.init, initargs=(tree,)) as ex:
            out.update(H=ex.submit(w.work, p).result(timeout=300))
    except Exception as e:
        out.update(pool_exception=type(e).__name__ + ": " + str(e)[:150])
print("JSON " + json.dumps(out), flush=True)
"""

WORKER = '''
import sys
def init(tree):
    sys.path.insert(0, tree)
def work(p):
    import camb
    return camb.get_results(p).hubble_parameter([0.5, 1.0]).tolist()
'''


def eft_keys():
    import harness

    keys = harness.parse_ini(FR_INI)
    return {k: v for k, v in keys.items() if k.startswith("EFT")}


def run(tree, mode, eft, ini=None):
    env = dict(os.environ, OMP_NUM_THREADS="2", PYTHONPATH=str(WD))
    p = subprocess.run([sys.executable, "-c", CHILD, tree, mode, json.dumps(eft)] + ([ini] if ini else []),
                       capture_output=True, text=True, env=env, cwd=WD, timeout=600)
    res = None
    for ln in p.stdout.splitlines():
        if ln.startswith("JSON "):
            res = json.loads(ln[5:])
    return p.returncode, res, p.stdout.strip().splitlines()[-1:] if p.stdout.strip() else [], p.stderr[-300:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=f"{ROOT}/trees/A0")
    ap.add_argument("--ini", default=None, help="build base params with read_ini of this file")
    ap.add_argument("--modes", default="state,ref,gr_pickle,deepcopy_results,set_params_cambdata,workaround,pool")
    a = ap.parse_args()
    WD.mkdir(parents=True, exist_ok=True)
    (WD / "_verify_pool_worker_attr.py").write_text(WORKER)
    eft = eft_keys()
    print("tree:", a.tree, "ini:", a.ini)
    print("EFT keys:", eft)
    for mode in a.modes.split(","):
        rc, res, last, err = run(a.tree, mode, eft, a.ini)
        st = "OK" if rc == 0 else (f"signal {-rc}" if rc < 0 else f"rc={rc}")
        print(f"{mode:20s} {st:10s} {json.dumps(res) if res else last} {'' if rc == 0 else err.strip()[-200:]}")


if __name__ == "__main__":
    main()

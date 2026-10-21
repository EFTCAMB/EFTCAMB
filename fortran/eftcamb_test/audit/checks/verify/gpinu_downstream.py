"""Verifier: downstream effect of the tripled gpinu_tot in the stability cache.

For each model, fills the stability cache with __eftcamb_stability_MOD_eftteststability on a log grid of a, in two
trees (code under audit and the fixed tree), and prints max relative change of the stability inputs.

Usage: python gpinu_downstream.py OUT.json TREE_BUGGY TREE_FIXED model.ini [model.ini ...]
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

FIELDS = ["adotoa", "Hdot", "Hdotdot", "Hdotdotdot", "EFTcdot", "EFTLambdadot", "EFTpiA1", "EFTpiA2", "EFTpiB1", "EFTpiB2", "EFTpiC",
          "EFTpiD1", "EFTpiD2", "EFT_kinetic", "EFT_gradient", "EFT_mu1", "EFT_mu2", "posbound1", "posbound2",
          "ghostpos", "tachpos"]

WORKER = r"""
import ctypes, json, sys
from ctypes import byref, c_double, c_int
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
from camb.baseconfig import camblib as lib
from camb.eftcamb import EFTCAMB_timestep_cache
cp = camb.read_ini(sys.argv[2], no_validate=True)
cp.EFTCAMB.EFTCAMB_skip_stability = True
res = camb.get_background(cp)
P = res.Params
test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
buf = (ctypes.c_char * (1 << 18))()
names = json.loads(sys.argv[3])
off = {n: getattr(EFTCAMB_timestep_cache, n).offset for n in names}
a = np.logspace(-8, 0, 161)
out = {n: [] for n in names}
for x in a:
    test(byref(c_double(x)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf)
    raw = bytes(buf[: max(off.values()) + 8])
    for n in names:
        out[n].append(float(np.frombuffer(raw, dtype=np.float64, count=1, offset=off[n])[0]))
out["a"] = a.tolist()
print("RESULT " + json.dumps(out))
"""


def run(tree, ini):
    with tempfile.TemporaryDirectory() as d:
        w = Path(d) / "m.ini"
        w.write_text(f"DEFAULT({ini})\noutput_root = {d}/out\nfeedback_level = 0\nEFTCAMB_write_background = F\n")
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        with harness.slot(1):
            p = subprocess.run([sys.executable, "-c", WORKER, tree, str(w), json.dumps(FIELDS)], capture_output=True,
                               text=True, cwd=d, env=env, preexec_fn=lambda: os.nice(10), timeout=900)
        res = [json.loads(line[7:]) for line in p.stdout.splitlines() if line.startswith("RESULT ")]
        if not res:
            return {"error": (p.stdout[-1500:] + p.stderr[-1500:])}
        return res[0]


def main():
    out_json, tb, tf = sys.argv[1:4]
    allres = {}
    for ini in sys.argv[4:]:
        rb, rf = run(tb, ini), run(tf, ini)
        name = Path(ini).stem
        allres[name] = {"buggy": rb, "fixed": rf}
        print(f"== {name}")
        if "error" in rb or "error" in rf:
            print("  ERROR", (rb.get("error") or rf.get("error"))[-400:])
            continue
        a = np.array(rb["a"])
        for n in FIELDS:
            b, f = np.array(rb[n]), np.array(rf[n])
            scale = np.maximum(np.abs(f), 1e-300)
            rel = np.abs(b - f) / scale
            rel[np.abs(f) == 0] = np.where(b[np.abs(f) == 0] == 0, 0.0, np.inf)
            i = int(np.nanargmax(rel)) if np.any(np.isfinite(rel)) else 0
            sign = int(np.sum(np.sign(b) != np.sign(f)))
            print(f"   {n:14s} max rel change {rel[i]:10.3e} at a={a[i]:8.1e}   sign flips over grid: {sign}")
    Path(out_json).write_text(json.dumps(allres))


if __name__ == "__main__":
    main()

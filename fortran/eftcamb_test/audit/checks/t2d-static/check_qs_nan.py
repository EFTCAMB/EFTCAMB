"""Scan every test model for NaN / inf in the quasi-static outputs (QS_mu, QS_sigma, QS_mu_gen, QS_sigma_gen) and,
for reference, a few other per-k output fields, using EFTCAMB.get_scale_evolution after camb.get_background.

Usage: python check_qs_nan.py OUT.json [TREE]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
PAR = REPO / "fortran/eftcamb_test/parameters"
SCR = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit")
FIELDS = ["QS_mu", "QS_sigma", "QS_mu_gen", "QS_sigma_gen", "mu", "gamma", "sigma_eft", "Psi", "Phi", "EFTISW",
          "EFTLensing", "delta_DE_eff"]
CHILD = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb, numpy as np
fields = json.loads(sys.argv[3])
p = camb.read_ini(sys.argv[2], no_validate=True)
r = camb.get_background(p)
q = np.array([1e-3, 0.01, 0.1, 1.0]); a = np.array([0.02, 0.1, 0.5, 1.0])
names, arr = p.EFTCAMB.get_scale_evolution(r, q, a)
out = {}
for k in fields:
    v = np.asarray(arr[k]).ravel()
    out[k] = int((~np.isfinite(v)).sum())
out["n"] = int(v.size)
print("JSON" + json.dumps(out))
"""


def run(tree, ini, tmp):
    w = tmp / f"{ini.stem}.ini"
    w.write_text(f"DEFAULT({ini})\nfeedback_level = 0\nEFTCAMB_feedback_level = 0\nEFTCAMB_write_background = F\n"
                 f"output_root = {tmp / ini.stem}\n")
    try:
        r = subprocess.run([sys.executable, "-c", CHILD, str(tree), str(w), json.dumps(FIELDS)], capture_output=True,
                           text=True, cwd=tmp, env=dict(os.environ, OMP_NUM_THREADS="2"), timeout=1200)
    except subprocess.TimeoutExpired:
        return {"status": "timeout"}
    js = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
    if not js:
        return {"status": f"failed rc={r.returncode}", "tail": (r.stdout + r.stderr)[-300:]}
    return dict(status="ok", **json.loads(js[-1][4:]))


def main():
    out = Path(sys.argv[1])
    tree = Path(sys.argv[2]) if len(sys.argv) > 2 else SCR / "trees/A0"
    inis = [i for i in sorted(PAR.glob("[2-9]*.ini"))]
    tmp = Path(tempfile.mkdtemp(prefix="qsnan_"))
    with ThreadPoolExecutor(4) as ex:
        res = dict(zip([i.stem for i in inis], ex.map(lambda i: run(tree, i, tmp), inis)))
    bad = 0
    for m, d in res.items():
        if d["status"] != "ok":
            print(f"{m:30s} {d['status']}")
            continue
        nonfin = {k: d[k] for k in FIELDS if d[k]}
        bad += bool(nonfin)
        print(f"{m:30s} non-finite (of {d['n']}): {nonfin if nonfin else 'none'}")
    print(f"models: {len(res)}, with non-finite outputs: {bad}")
    out.write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

"""Repro: re-using the parameters held by a CAMBdata (results.Params, or a copy of it) aborts the process.

CAMBdata_SetParams copies P into this%CP (so an already-filled eft_par_cache comes along) and then does
    allocate(this%CP%eft_par_cache%grhormass(max_nu), this%CP%eft_par_cache%nu_masses(max_nu))
without checking allocation (fortran/results.f90:580) -> "Attempting to allocate already allocated variable".
CAMB v2's results.get_BAO(z) defaults to self.Params.copy(), so the standard API call crashes for every EFT model.

Usage: python eft-par-cache-realloc-crash.py [--tree PATH]
Prints, for each scenario, observed (exit code / message) vs expected (same result as GR or as a fresh run).
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
PARAMS = Path(__file__).resolve().parents[2] / "parameters"

CHILD = r'''
import sys
sys.path.insert(0, sys.argv[1])
import camb
mode, ini, eft = sys.argv[2], sys.argv[3], sys.argv[4] == "1"
if eft:
    cp = camb.read_ini(ini, no_validate=True)
else:
    cp = camb.set_params(H0=70, ombh2=0.0226, omch2=0.112, mnu=0.06)
res = camb.get_background(cp)
if mode == "get_BAO":
    print("RESULT", res.get_BAO([0.5, 1.0])[:, 1].tolist(), flush=True)
elif mode == "get_background(results.Params)":
    print("RESULT", camb.get_background(res.Params).hubble_parameter([0.5, 1.0]).tolist(), flush=True)
elif mode == "get_background(results.Params.copy())":
    print("RESULT", camb.get_background(res.Params.copy()).hubble_parameter([0.5, 1.0]).tolist(), flush=True)
elif mode == "reference":
    print("RESULT", res.hubble_parameter([0.5, 1.0]).tolist(), flush=True)
'''


def run(tree, mode, ini, eft):
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", CHILD, tree, mode, ini, "1" if eft else "0"], capture_output=True,
                       text=True, env=env, cwd=os.path.dirname(ini))
    res = [ln for ln in p.stdout.splitlines() if ln.startswith("RESULT")]
    err = [ln for ln in (p.stdout + p.stderr).splitlines() if "Fortran runtime error" in ln or "At line" in ln]
    return p.returncode, (res[0][7:] if res else None), err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=f"{ROOT}/trees/A0")
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix="realloc_")
    ini = os.path.join(tmp, "fR.ini")
    Path(ini).write_text(f"DEFAULT({PARAMS / '4_Des_F_R_LCDM_1.ini'})\noutput_root = {tmp}/fR\n"
                         "EFTCAMB_write_background = F\nfeedback_level = 0\n")
    print(f"tree: {a.tree}")
    rc, ref, _ = run(a.tree, "reference", ini, True)
    print(f"designer f(R) fresh run, H(z=0.5,1) = {ref} (rc={rc})")
    for eft in (False, True):
        for mode in ("get_BAO", "get_background(results.Params)", "get_background(results.Params.copy())"):
            rc, out, err = run(a.tree, mode, ini, eft)
            label = "EFTflag=3 designer f(R)" if eft else "GR (EFTflag=0)       "
            status = "OK" if rc == 0 else f"CRASH rc={rc}"
            print(f"{label} {mode:40s} observed: {status} {out or ' '.join(err)}")
            if eft:
                print(f"{'':22s} {'':40s} expected: rc=0 and a result (as for GR / fresh run {ref})")


if __name__ == "__main__":
    main()

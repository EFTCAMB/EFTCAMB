"""Verify (mass-stability-verdict-unreliable): real-code verdicts with EFT_mass_stability forced on, read from the
Fortran feedback (EFTCAMB_feedback_level = 2 prints 'Mass instability. mu_i/H^2 = ...' at the first failing time).

Run: python massverd_verdicts.py tree model [model ...]
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

PAR = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"
WORKER = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
try:
    camb.get_background(cp)
    print("VERDICT stable", flush=True)
except Exception as e:
    print("VERDICT unstable", repr(e)[:120], flush=True)
"""


def verdict(tree, model, mass):
    with tempfile.TemporaryDirectory() as d:
        ini = Path(d) / "m.ini"
        ini.write_text(f"DEFAULT({PAR}/{model}.ini)\nfeedback_level = 2\n"
                       f"EFT_mass_stability = {'T' if mass else 'F'}\nEFTCAMB_write_background = F\n"
                       f"output_root = {d}/out\n")
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        with harness.slot(1):
            p = subprocess.run([sys.executable, "-c", WORKER, tree, str(ini)], capture_output=True, text=True,
                               cwd=d, env=env, timeout=600, preexec_fn=lambda: os.nice(10))
    out = p.stdout + p.stderr
    v = [ln for ln in out.splitlines() if ln.startswith("VERDICT")]
    mass_lines = [ln.strip() for ln in out.splitlines() if "Mass instability" in ln]
    other = [ln.strip() for ln in out.splitlines() if "Physical instability" in ln and "Mass" not in ln]
    nan = [ln.strip() for ln in out.splitlines() if "Nan" in ln or "NaN" in ln]
    return {"verdict": v[0].split()[1] if v else "?", "mass_msgs": mass_lines[:2], "other_msgs": other[:2],
            "nan_msgs": nan[:3]}


if __name__ == "__main__":
    tree = sys.argv[1]
    res = {}
    for m in sys.argv[2:]:
        r_off = verdict(tree, m, False)
        r_on = verdict(tree, m, True)
        res[m] = {"mass_off": r_off, "mass_on": r_on}
        print(f"{m:28s} mass=F: {r_off['verdict']:9s} mass=T: {r_on['verdict']:9s} "
              f"{(r_on['mass_msgs'] or r_on['nan_msgs'] or r_on['other_msgs'] or [''])[0]}")
    print("RESULT " + json.dumps(res))

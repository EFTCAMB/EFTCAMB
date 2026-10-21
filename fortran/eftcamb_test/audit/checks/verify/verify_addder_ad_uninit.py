"""Debug-build (AD: -finit-real=snan -ffpe-trap=invalid) check of the uninitialised gamma2PP/gamma3PP/gamma3PPP dummies
in EFTCAMBFullModelComputeAdditionalDer: run a full-mapping model with EFT_mass_stability = T vs F through
camb.get_background (stability sampler only). Usage: python verify_addder_ad_uninit.py TREE INI OUTDIR
"""

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness  # noqa: E402

CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.EFTCAMB.EFT_mass_stability = sys.argv[3] == "T"
pars.EFTCAMB.EFTCAMB_feedback_level = int(sys.argv[4])
try:
    camb.get_background(pars)
    print("RESULT ok")
except Exception as e:
    print("RESULT exception", type(e).__name__, str(e).strip().splitlines()[-1])
"""

tree, ini, outdir = sys.argv[1], sys.argv[2], Path(sys.argv[3])
outdir.mkdir(parents=True, exist_ok=True)
w = outdir / f"w_{Path(ini).stem}.ini"
w.write_text(f"DEFAULT({ini})\nfeedback_level = 0\nEFTCAMB_write_background = F\noutput_root = {outdir / 'x'}\n")
for flag in ("F", "T"):
    with harness.slot(1):
        r = subprocess.run([sys.executable, "-c", CHILD, tree, str(w), flag, os.environ.get("FB", "0")], capture_output=True, text=True,
                           cwd=outdir, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900,
                           preexec_fn=lambda: os.nice(10))
    res = [ln for ln in r.stdout.splitlines() if ln.startswith("RESULT")]
    print(f"EFT_mass_stability={flag}: rc={r.returncode} {res[-1] if res else ''}")
    if r.returncode != 0 or not res or os.environ.get("FB"):
        print("   stdout tail:", "\n   ".join([x for x in r.stdout.strip().splitlines() if "nstab" in x or "NaN" in x or "nan" in x][-8:]))
        print("   stderr tail:", "\n   ".join(r.stderr.strip().splitlines()[-14:]))

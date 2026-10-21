"""Verifier: EFTCAMB get_scale_evolution after transfer-only run. Independent path: set params via
set_params/python API (no read_ini), vary model and call orders; compare with upstream get_time_evolution."""
import signal
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "getevol"
WORKER = r"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, sys.argv[1])
import camb, numpy as np
mode, call, ini = sys.argv[2], sys.argv[3], sys.argv[4]
if ini == "none":
    cp = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965)
else:
    cp = camb.read_ini(ini, no_validate=True)
    cp.SourceWindows = []
    cp.WantTensors = False
    cp.DoLensing = False
    cp.InitPower.r = 0.0
if mode == "transfer":
    cp.WantCls = False; cp.WantScalars = False
cp.WantTransfer = True
cp.set_matter_power(redshifts=[0.0], kmax=1.0)
r = camb.get_results(cp)
print("WantScalars after run:", r.Params.WantScalars, flush=True)
if call == "upstream":
    ev = r.get_time_evolution([0.05], np.array([100., 1000., 10000.]), ["delta_cdm"])
    print("upstream ok", np.asarray(ev).ravel().tolist(), flush=True)
else:
    names, arr = cp.EFTCAMB.get_scale_evolution(r, [0.05], np.array([0.1, 0.5, 1.0]))
    key = "pi" if "pi" in (arr.dtype.names or names) else names[0]
    print("eft ok", key, np.asarray(arr[key])[0].tolist(), flush=True)
"""

def case(tree, mode, call, ini):
    OUT.mkdir(parents=True, exist_ok=True)
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-I", "-c", WORKER, str(harness.tree_path(tree)), mode, call, ini],
                           capture_output=True, text=True, timeout=600, cwd=OUT)
    rc = p.returncode
    sig = f" ({signal.Signals(-rc).name})" if rc < 0 else ""
    return f"rc={rc}{sig} " + " | ".join(ln for ln in p.stdout.splitlines()[-2:]) + (" ERR:" + p.stderr.strip().splitlines()[-1] if rc and p.stderr.strip() else "")

ini_rph = str(harness.PARAMS_DIR / "3_RPH_MK_const_1.ini")
ini_fr = str(harness.PARAMS_DIR / sorted(p.name for p in harness.PARAMS_DIR.glob("*fR*.ini"))[0]) if list(harness.PARAMS_DIR.glob("*fR*.ini")) else ini_rph
cases = [("A0", "transfer", "eft", ini_rph), ("A0", "transfer", "eft", ini_fr), ("A0", "transfer", "eft", "none"),
         ("A0", "cls", "eft", "none"), ("A0", "transfer", "upstream", ini_rph), ("A0", "transfer", "upstream", "none"),
         ("B0", "transfer", "upstream", "none"), ("B0", "cls", "upstream", "none"), ("A0", "cls", "upstream", "none")]
for c in cases:
    print(c[0], c[1], c[2], Path(c[3]).name, "->", case(*c), flush=True)

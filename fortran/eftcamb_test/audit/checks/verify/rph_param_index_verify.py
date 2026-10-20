"""Independent check of rph-param-index-tensor-ode: several ODE-block configurations, A0 and AD (snan init) trees.

For each case prints num_params and every (name, latex, value) from the Python bridge, plus the expected list
built from the ini (block order wDE, M, K, B, T, M_ODE, K_ODE, B_ODE, T_ODE).
"""

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

REPO = Path(__file__).resolve().parents[5]
S = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
         "scratchpad/audit")
OUT = S / "work/verify/rphidx"
INI = REPO / "fortran/eftcamb_test/parameters/3_RPH_MK_power_1.ini"
CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
e = pars.EFTCAMB
n = e.num_params()
print("num_params", n)
vals = list(e.param_values())
print("values", vals)
if sys.argv[3] == "1":
    print("names", list(e.param_names()))
"""
CASES = {
    "massPODE_only": "RPHmassPmodel_ODE = 1\nRPHmassP_ODE0 = 0.123\n",
    "tensorODE_only": "RPHtensormodel_ODE = 1\nRPHtensor_ODE0 = 0.456\n",
    "both_ODE": "RPHmassPmodel_ODE = 1\nRPHmassP_ODE0 = 0.123\nRPHtensormodel_ODE = 1\nRPHtensor_ODE0 = 0.456\n",
    "massPODE_kinODE": "RPHmassPmodel_ODE = 1\nRPHmassP_ODE0 = 0.123\nRPHkineticitymodel_ODE = 1\n"
                       "RPHkineticity_ODE0 = 0.789\n",
}
for tree in ("A0", "AD"):
    for name, extra in CASES.items():
        w = OUT / f"{tree}_{name}.ini"
        w.write_text(f"DEFAULT({INI})\nfeedback_level = 0\n{extra}EFTCAMB_write_background = F\n"
                     f"output_root = {OUT / 'out'}\n")
        with harness.slot(1):
            r = subprocess.run([sys.executable, "-c", CHILD, str(S / "trees" / tree), str(w),
                                "1" if tree == "A0" else "0"], capture_output=True, text=True, cwd=OUT,
                               env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=600)
        print(f"== {tree} {name} rc={r.returncode}")
        print(r.stdout.strip())
        if r.returncode:
            print(r.stderr[-1500:])

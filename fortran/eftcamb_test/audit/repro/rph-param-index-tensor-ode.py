"""Repro: RPH parameter getters use RPH_Tensor_ODE instead of RPH_PlanckMass_ODE for the M_ODE block size.

Code: fortran/eftcamb/07f_designer_models/007p2_RPH.f90:610 (ParameterNames), :731 (ParameterNamesLatex),
:852 (ParameterValues), with RPHusealphaM = F:
        NMo = NT + self%RPH_Tensor_ODE%parameter_number          ! should be RPH_PlanckMass_ODE
while ComputeParametersNumber (:515) and InitModelParameters (:413) use RPH_PlanckMass_ODE.  When the two ODE
functions have different parameter counts every parameter after NT is mapped to the wrong block or to none
(the intent(out) name/value is then left undefined).

Case: 3_RPH_MK_power_1.ini + RPHmassPmodel_ODE = 1 (constant, 1 parameter RPHmassP_ODE0 = 0.123), Tensor_ODE = 0.
Expected: param_names()[-1] == 'RPHmassP_ODE0' and param_values()[-1] == 0.123.
Usage: /Users/marco/Software/miniforge3/bin/python rph-param-index-tensor-ode.py [TREE]
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0"
)
INI = REPO / "fortran/eftcamb_test/parameters/3_RPH_MK_power_1.ini"
CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
e = pars.EFTCAMB
print("num_params", e.num_params())
print("names ", list(e.param_names()))
print("values", list(e.param_values()))
"""


def main():
    tmp = Path(tempfile.mkdtemp(prefix="rphidx_"))
    w = tmp / "w.ini"
    w.write_text(f"DEFAULT({INI})\nfeedback_level = 0\nRPHmassPmodel_ODE = 1\nRPHmassP_ODE0 = 0.123\n"
                 f"EFTCAMB_write_background = F\noutput_root = {tmp / 'out'}\n")
    r = subprocess.run([sys.executable, "-c", CHILD, str(TREE), str(w)], capture_output=True, text=True, cwd=tmp,
                       env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=600)
    print(f"tree: {TREE}  rc={r.returncode}")
    print(r.stdout.strip()[-3000:])
    if r.returncode:
        print(r.stderr[-2000:])
    print("expected: last name 'RPHmassP_ODE0', last value 0.123")


if __name__ == "__main__":
    main()

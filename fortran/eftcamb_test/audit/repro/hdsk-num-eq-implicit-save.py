"""Repro: Horndeski background solver keeps num_eq=6 across runs (implicit SAVE of `integer :: num_eq = 5`).

fortran/eftcamb/08f_full_models/008p0_Horndeski.f90:2381  `integer :: num_eq = 5` is initialised in its declaration,
so it has the SAVE attribute; line 2426 sets it to 6 when Horndeski_evolve_hubble (default T) or designer_V is used and
nothing ever resets it.  A later Horndeski run with Horndeski_evolve_hubble = F in the same process then integrates a
6-component system whose y(6) and ydot(6) are never set.

Procedure (one Python process, tree A0):
  run 1: 5_hdsk_TG_1.ini (evolve_hubble = F)               -> reference
  run 2: same ini with Horndeski_evolve_hubble = T          -> sets the saved num_eq to 6
  run 3: 5_hdsk_TG_1.ini again (evolve_hubble = F)          -> should be bit-identical to run 1
Expected: run 3 == run 1 (max rel. difference 0). Observed: printed below.
A second fresh process runs run1 twice to show that a plain repeat is reproducible.

Usage: /Users/marco/Software/miniforge3/bin/python hdsk-num-eq-implicit-save.py [TREE]
"""

from __future__ import annotations

import json
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
INI = REPO / "fortran/eftcamb_test/parameters/5_hdsk_TG_1.ini"

CHILD = r"""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
seq = json.loads(sys.argv[2])
z = np.concatenate([np.linspace(0, 3, 31), np.logspace(0.7, 3.5, 30)])
a = np.logspace(-6, 0, 61)
out = []
for ini in seq:
    pars = camb.read_ini(ini, no_validate=True)
    # note: setting WantCls = WantTransfer = False here makes get_background segfault (also in upstream B0 with this
    # ini), so keep the ini flags; get_background only computes the background either way.
    res = camb.get_background(pars)
    H = np.asarray(res.hubble_parameter(z))
    fields, arr = pars.EFTCAMB.get_eft_functions(res, a)
    out.append(dict(H=H.tolist(), adotoa=arr["adotoa"].tolist(), Hdot=arr["Hdot"].tolist(),
                    OmegaV=arr["EFTOmegaV"].tolist(), c=arr["EFTc"].tolist(), tau0=float(res.tau0)))
print("JSON" + json.dumps(out))
"""


def wrapper(tmp: Path, name: str, extra: dict[str, str]) -> str:
    p = tmp / f"{name}.ini"
    lines = [f"DEFAULT({INI})", "feedback_level = 0", "EFTCAMB_write_background = F",
             f"output_root = {tmp / name}"]
    lines += [f"{k} = {v}" for k, v in extra.items()]
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def run(seq: list[str], tmp: Path) -> list[dict]:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    r = subprocess.run([sys.executable, "-c", CHILD, str(TREE), json.dumps(seq)], capture_output=True, text=True,
                       env=env, cwd=tmp, timeout=900)
    line = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
    if not line:
        print("child failed rc=", r.returncode, r.stdout[-2000:], r.stderr[-2000:])
        sys.exit(1)
    return json.loads(line[-1][4:])


def maxrel(x, y):
    import numpy as np

    x, y = np.asarray(x), np.asarray(y)
    return float(np.max(np.abs(x - y) / np.maximum(np.abs(x), 1e-300)))


def main():
    tmp = Path(tempfile.mkdtemp(prefix="hdsk_numeq_"))
    base = wrapper(tmp, "tg1_F", {"Horndeski_evolve_hubble": "F"})
    evol = wrapper(tmp, "tg1_T", {"Horndeski_evolve_hubble": "T"})
    seq = run([base, evol, base], tmp)
    ctrl = run([base, base], tmp)
    print(f"tree: {TREE}")
    print("control (fresh process, F then F): max rel diff  H(z) = %.3e  adotoa = %.3e  Hdot = %.3e  tau0 = %.3e"
          % (maxrel(ctrl[0]["H"], ctrl[1]["H"]), maxrel(ctrl[0]["adotoa"], ctrl[1]["adotoa"]),
             maxrel(ctrl[0]["Hdot"], ctrl[1]["Hdot"]), abs(ctrl[0]["tau0"] / ctrl[1]["tau0"] - 1)))
    print("sequence F, T, F: run3 vs run1 max rel diff  H(z) = %.3e  adotoa = %.3e  Hdot = %.3e  OmegaV = %.3e  "
          "EFTc = %.3e  tau0 = %.3e"
          % (maxrel(seq[0]["H"], seq[2]["H"]), maxrel(seq[0]["adotoa"], seq[2]["adotoa"]),
             maxrel(seq[0]["Hdot"], seq[2]["Hdot"]), maxrel(seq[0]["OmegaV"], seq[2]["OmegaV"]),
             maxrel(seq[0]["c"], seq[2]["c"]), abs(seq[0]["tau0"] / seq[2]["tau0"] - 1)))
    print("expected: run3 identical to run1 (0.000e+00), as in the control")


if __name__ == "__main__":
    main()

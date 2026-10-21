"""End-to-end ghost-verdict scan (verifier for finding eft-kinetic-cross-term-extra-h0).

Uses EFTCAMB's own stability decision through camb.get_background (no spectra, avoids the get_eft_functions path).
RPH with alpha_M = m Omega_DE, alpha_B = b Omega_DE, alpha_K = k Omega_DE^2, alpha_T = 0. All of D scales as Omega_DE^2:
    correct no-ghost   D_true = k + 6 b^2
    code as written    D_code = k + 6 b^2 + 3 m^2 - 6 m b   (h0 -> 0 limit of the double-h0 bracket)
so the code accepts ghosts when D_true < 0 < D_code (b = 0, -3 m^2 < k < 0) and rejects healthy theories when
D_code < 0 < D_true (b > m/2).  The fixed tree (t6 fixkin) must follow D_true.
Each (tree, point) runs in its own subprocess under one harness CPU slot; feedback_level 3 makes EFTCAMB print the reason.
"""

import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ROOT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45"
            "/scratchpad/audit")
OUT = ROOT / "work/verify/ekin"
PAR = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")
TREES = {"A0": ROOT / "trees/A0", "fixkin": ROOT / "work/t6-stability/fixkin"}

# (m, b, k)
POINTS = [(0.5, 0.0, 0.1), (0.5, 0.0, -0.2), (0.5, 0.0, -0.5), (0.5, 0.0, -0.73), (0.5, 0.0, -0.77),
          (0.5, 0.0, -1.0), (0.3, 0.0, -0.25), (0.3, 0.0, -0.29), (0.0, 0.0, -0.05),
          (0.5, 0.5, -0.5), (0.5, 0.5, -1.0), (0.5, 0.5, -2.0), (1.0, 1.0, -4.0), (1.0, 1.0, -2.0)]
TAG = ""
if len(sys.argv) > 1 and sys.argv[1] == "braid":
    # false-rejection search: 3 m^2 - 6 m b < 0, k at the midpoint between D_code = 0 and D_true = 0
    TAG = "_braid"
    POINTS = []
    for m in (0.1, 0.3, -0.1, -0.3):
        for f in (0.75, 1.0, 2.0, -1.0):
            b = f * m
            if 3 * m**2 - 6 * m * b < 0:
                POINTS.append((m, b, round(-6 * b**2 - (3 * m**2 - 6 * m * b) / 2, 4)))

WORKER = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
try:
    camb.get_background(cp)
    print("VERDICT stable", flush=True)
except Exception as e:
    print("VERDICT unstable", repr(e)[:200], flush=True)
"""


def ini_text(m: float, b: float, k: float) -> str:
    s = (f"DEFAULT({PAR}/3_RPH_MK_OmegaDE_1.ini)\nfeedback_level = 3\n"
         "EFTCAMB_write_background = F\nget_scalar_cls = F\nget_tensor_cls = F\nget_transfer = F\ndo_lensing = F\n"
         f"RPHalphaMmodel_ODE = {2 if m != 0 else 0}\nRPHalphaM_ODE0 = {m}\n"
         f"RPHkineticitymodel_ODE = 3\nRPHkineticity_ODE0 = {k}\nRPHkineticity_ODEExp = 2.0\n")
    if b != 0:
        s += f"RPHbraidingmodel_ODE = 2\nRPHbraiding_ODE0 = {b}\n"
    return s


def run(tree: Path, m: float, b: float, k: float) -> dict:
    with tempfile.TemporaryDirectory() as d:
        ini = Path(d) / "m.ini"
        ini.write_text(ini_text(m, b, k))
        env = dict(os.environ, OMP_NUM_THREADS="1", GFORTRAN_UNBUFFERED_ALL="1")
        env.pop("PYTHONPATH", None)
        with harness.slot(1):
            p = subprocess.run([sys.executable, "-c", WORKER, str(tree), str(ini)], capture_output=True, text=True,
                               cwd=d, env=env, timeout=600, preexec_fn=lambda: os.nice(10))
    out = p.stdout + p.stderr
    v = re.findall(r"VERDICT (\w+)", out)
    reasons = sorted(set(re.findall(r"(?:Physical|Mathematical) instability: ([^\n=]+)", out)))
    kin = re.findall(r"ghost instability\. Kinetic term =\s*(\S+)", out)
    return {"verdict": v[0] if v else f"rc={p.returncode}", "reasons": [r.strip() for r in reasons],
            "first_ghost_kinetic": kin[0] if kin else None, "tail": out[-300:] if not v else ""}


rows = []
for m, b, k in POINTS:
    d_true = k + 6 * b**2
    d_code = k + 6 * b**2 + 3 * m**2 - 6 * m * b
    row = {"m": m, "b": b, "k": k, "D_true": d_true, "D_code": d_code}
    for name, tree in TREES.items():
        row[name] = run(tree, m, b, k)
    rows.append(row)
    print(f"m={m:.2f} b={b:.2f} k={k:+.2f}  D_true={d_true:+.3f} D_code={d_code:+.3f} | "
          + " | ".join(f"{n}: {row[n]['verdict']:8s} {row[n]['reasons']}" for n in TREES), flush=True)
(OUT / f"verdict_scan{TAG}.json").write_text(json.dumps(rows, indent=1))

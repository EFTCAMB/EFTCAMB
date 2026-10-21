"""Track 5: f(R) designer (B0 = 1e-4) turn-on / RGR-threshold convergence and scalaron oscillation diagnostics.

Runs sequentially (one CAMB process at a time):
  grt1e-10, grt1e-11 : EFTCAMB_GR_threshold = 1e-10, 1e-11 (default 1e-8, batch covers 1e-7 and 1e-9)
  evol_AB2           : accuracy_boost = 2 with get_scale_evolution at k = 0.1, 1, 2
  evol_pidfix        : tree with the EFTpiDfunction fix, same evolution output
Outputs under $WORK/fr/.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "t5-robust"
INI = harness.PARAMS_DIR / "4_Des_F_R_LCDM_2.ini"
JOBS = [
    ("grt1e-10", "A0", {"EFTCAMB_GR_threshold": "1e-10"}, None),
    ("grt1e-11", "A0", {"EFTCAMB_GR_threshold": "1e-11"}, None),
    ("evol_AB2", "A0", {"accuracy_boost": "2"}, [0.1, 1.0, 2.0]),
    ("evol_pidfix", W / "A0_pidfix", {}, [0.1, 1.0, 2.0]),
]
only = sys.argv[1].split(",") if len(sys.argv) > 1 else None
for name, tree, ov, ks in JOBS:
    if only and name not in only:
        continue
    out = W / "fr" / f"{name}.npz"
    if out.with_suffix(".json").exists():
        continue
    st = harness.run_model(tree, INI, out, overrides=ov, evol_k=ks, threads=2)
    print(name, st["status"], st["elapsed"], flush=True)

# python-attribute toggle of the Horndeski evolution flag (model keeps the flag it read from the ini)
if not only or "pyset_dphi" in only:
    out = W / "fr" / "TG2_pyset_dphiT.npz"
    if not out.with_suffix(".json").exists():
        st = harness.run_model("A0", harness.PARAMS_DIR / "5_hdsk_TG_2.ini", out,
                               python_set={"EFTCAMB.EFTCAMB_evolve_delta_phi": True}, threads=2)
        print("TG2_pyset_dphiT", st["status"], st["elapsed"], flush=True)

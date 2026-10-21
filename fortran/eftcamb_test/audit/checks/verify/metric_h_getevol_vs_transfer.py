"""Cross-check get_scale_evolution against the transfer functions for EFTCAMB_evolve_metric_h = F / T.

The k-mode trajectories from get_scale_evolution (T run) differ from the F run by O(1) factors at late times while the
transfer-level P(k) agrees to 1e-5; this script checks, in one process per flag, (1) delta_c(a=1) from the transfer
functions vs get_scale_evolution, (2) whether get_scale_evolution depends on the number of requested output times.

Usage: python metric_h_getevol_vs_transfer.py [tree]   (spawns one subprocess per flag)
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "metric_h"
KS = [1e-3, 1e-2, 0.1, 1.0]


def child(tree: str, flag: str) -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    import numpy as np

    sys.path.insert(0, str(harness.tree_path(tree)))
    import camb

    ini = OUT / f"gv_{tree}_{flag}.ini"
    ini.write_text(
        f"DEFAULT({harness.PARAMS_DIR / '5_hdsk_TG_1.ini'})\noutput_root = {OUT / ('gv_' + flag)}\n"
        f"EFTCAMB_write_background = F\nEFTCAMB_evolve_metric_h = {flag}\n"
    )
    cp = camb.read_ini(str(ini), no_validate=True)
    if cp.EFTCAMB.EFTflag == 0:
        cp.EFTCAMB.initialize_parameters(cp, harness.parse_ini(ini), print_header=False)
    res = camb.get_results(cp)
    tr = res.get_matter_transfer_data()
    kh = tr.transfer_data[0, :, -1] * cp.H0 / 100  # k in 1/Mpc
    dc = tr.transfer_data[camb.model.Transfer_cdm - 1, :, -1]  # delta_c/k^2 at z = 0 (last redshift)
    out = {"evolve_metric_h": bool(cp.EFTCAMB.EFTCAMB_evolve_metric_h), "pk_redshifts": list(map(float, cp.Transfer.PK_redshifts[:cp.Transfer.PK_num_redshifts]))}
    out["transfer_dc_over_k2"] = [float(np.exp(np.interp(np.log(k), np.log(kh), np.log(np.abs(dc))))) for k in KS]
    for n in (300, 3, 1):
        a = np.logspace(-6, 0, n) if n > 1 else np.array([1.0])
        _, arr = cp.EFTCAMB.get_scale_evolution(res, KS, a)
        out[f"getevol_n{n}_clxc_a1"] = [float(arr["clxc"][i, -1]) for i in range(len(KS))]
        out[f"getevol_n{n}_z_a1"] = [float(arr["z"][i, -1]) for i in range(len(KS))]
        out[f"getevol_n{n}_etak_a1"] = [float(arr["etak"][i, -1]) for i in range(len(KS))]
    print("RESULT " + json.dumps(out))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--child":
        with harness.slot(1):
            child(sys.argv[2], sys.argv[3])
        sys.exit(0)
    tree = sys.argv[1] if len(sys.argv) > 1 else "A0"
    r = {}
    for flag in ("F", "T"):
        p = subprocess.run(
            [sys.executable, __file__, "--child", tree, flag], capture_output=True, text=True, timeout=1500,
            env=dict(os.environ, OMP_NUM_THREADS="1"),
        )
        line = [x for x in p.stdout.splitlines() if x.startswith("RESULT ")]
        if not line:
            raise SystemExit(p.stdout[-2000:] + p.stderr[-2000:])
        r[flag] = json.loads(line[0][7:])
    print(f"tree={tree}, k = {KS}")
    for key in r["F"]:
        if isinstance(r["F"][key], list):
            rat = [t / f if f else float("nan") for t, f in zip(r["T"][key], r["F"][key])]
            print(f"{key:28s} F: " + " ".join(f"{x:+.4e}" for x in r["F"][key]))
            print(f"{'':28s} T: " + " ".join(f"{x:+.4e}" for x in r["T"][key]) + "   T/F: " +
                  " ".join(f"{x:.4f}" for x in rat))
    (OUT / f"getevol_vs_transfer_{tree}.json").write_text(json.dumps(r, indent=1))

"""Independent check of t3-pi-ic-singular-denominator (adversarial verify).

Instead of the P(k) transfer grid, evolve single k modes with EFTCAMB's own per-mode evolution output
(get_scale_evolution -> EFTCAMB_quantities_evolution) on a fine k grid around k* = sqrt(-C/D1) at the turn-on,
and compare clxc(a=1) of RPH (alpha_K=0.1, alpha_M=m) against the same model with alpha_M=0, for:
  IC1 : EFT_IC_type=1 (code under audit, pi = -h0 E/(C + k^2 D1 + k^4 D2))
  IC2 : EFT_IC_type=2 (radiation-era series, finite for all k; type 4 -> `case default` is rejected at read time)
  IC3 : EFT_IC_type=3
A pole in the IC shows up as |dP/P| ~ 1/|k-k*| with a sign flip across k*, located where the denominator evaluated at
the turn-on time vanishes; a physical feature would be smooth in k and independent of the IC.

Usage: python piic_pole_scan.py <alpha_M> <outjson> [turn_on]
"""

import json
import os
import sys
import tempfile
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

import numpy as np  # noqa: E402

TREE = harness.TREES["A0"]
sys.path.insert(0, str(TREE))
import camb  # noqa: E402

assert Path(camb.__file__).resolve().is_relative_to(TREE.resolve()), camb.__file__
INI = harness.PARAMS_DIR / "3_RPH_MK_const_1.ini"
BASE = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "100", "num_redshiftwindows": "0",
        "get_transfer": "F"}


def params(over: dict) -> camb.CAMBparams:
    lines = [f"DEFAULT({INI})", "EFTCAMB_write_background = F", "output_root = /dev/null/x"]
    lines += [f"{k} = {v}" for k, v in {**BASE, **over}.items()]
    with tempfile.NamedTemporaryFile("w", suffix=".ini", delete=False) as f:
        f.write("\n".join(lines) + "\n")
    try:
        return camb.read_ini(f.name, no_validate=True)
    finally:
        os.unlink(f.name)


def evolve(over: dict, ks: np.ndarray, a: np.ndarray) -> dict:
    cp = params(over)
    res = camb.get_results(cp)
    _, arr = cp.EFTCAMB.get_scale_evolution(res, ks, a)
    return {f: np.asarray(arr[f]) for f in ("clxc", "pi", "pidot", "EFTpiC", "EFTpiD1", "EFTpiD2", "EFTpiE",
                                            "EFTpiA1", "EFTpiA2")}


def main() -> None:
    m = sys.argv[1]
    out = Path(sys.argv[2])
    ton = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01
    a = np.array([ton * (1 + 1e-9), ton * 1.0001, ton * 1.01, 0.02, 0.05, 0.1, 0.3, 1.0])
    with harness.slot(1):
        # 1. locate k* from the coefficients at the turn-on (coarse k grid; C, D1, D2 do not depend on k)
        probe = evolve({"RPHmassP0": m, "EFTCAMB_turn_on_time": ton}, np.array([0.01, 0.1, 1.0]), a)
        j0 = int(np.argmax(probe["EFTpiD1"][0] != 0))
        C, D1, D2 = (float(probe[f][0, j0]) for f in ("EFTpiC", "EFTpiD1", "EFTpiD2"))
        print(f"alpha_M={m} turn_on={ton}: at a={a[j0]:.8f} C={C:+.4e} D1={D1:+.4e} D2={D2:+.4e} "
              f"A1={probe['EFTpiA1'][0, j0]:+.3e} A2={probe['EFTpiA2'][0, j0]:+.3e}")
        if D2 != 0:
            disc = D1 * D1 - 4 * D2 * C
            k2s = [(-D1 + s * np.sqrt(disc)) / (2 * D2) for s in (1, -1)] if disc >= 0 else []
            k2s = [x for x in k2s if x > 0]
        else:
            k2s = [-C / D1] if C * D1 < 0 else []
        kstar = float(np.sqrt(k2s[0])) if k2s else float("nan")
        print(f"  k* = sqrt(-C/D1) = {kstar:.6g} 1/Mpc")
        if not np.isfinite(kstar):
            out.write_text(json.dumps({"alpha_M": m, "kstar": None}))
            return
        rel = np.concatenate([-np.logspace(-1, -5, 9), np.logspace(-5, -1, 9)])
        ks = np.sort(np.concatenate([kstar * (1 + rel), kstar * np.array([0.5, 0.7, 1.4, 2.0])]))
        runs = {}
        for name, ov in (("IC1", {}), ("IC2", {"EFT_IC_type": "2"}), ("IC3", {"EFT_IC_type": "3"}),
                         ("ref", {"RPHmassP0": "0.0"})):
            try:
                runs[name] = evolve({"RPHmassP0": m, "EFTCAMB_turn_on_time": ton, **ov}, ks, a)
            except Exception as e:  # noqa: BLE001
                print(f"  {name}: failed: {str(e)[:90]}")
                runs[name] = {f: np.full((len(ks), len(a)), np.nan) for f in ("clxc", "pi")}
    ref = runs["ref"]["clxc"][:, -1]
    rows = []
    print(f"  {'k/k*-1':>10s} {'dP/P IC1':>11s} {'dP/P IC2':>11s} {'dP/P IC3':>11s} {'pi_ini IC1':>11s} "
          f"{'pi_ini IC2':>11s} {'pi(a=1) IC1':>12s} {'pi(a=1) IC2':>12s} {'den(k) at ton':>13s}")
    for i, k in enumerate(ks):
        r1 = (runs["IC1"]["clxc"][i, -1] / ref[i]) ** 2 - 1
        r2 = (runs["IC2"]["clxc"][i, -1] / ref[i]) ** 2 - 1
        r3 = (runs["IC3"]["clxc"][i, -1] / ref[i]) ** 2 - 1
        den = C + k * k * D1 + k**4 * D2
        rows.append({"k": k, "k_over_kstar_m1": k / kstar - 1, "dP_IC1": r1, "dP_IC2": r2, "dP_IC3": r3,
                     "pi_ini_IC1": runs["IC1"]["pi"][i, j0], "pi_end_IC1": runs["IC1"]["pi"][i, -1],
                     "pi_ini_IC2": runs["IC2"]["pi"][i, j0], "pi_end_IC2": runs["IC2"]["pi"][i, -1], "den": den})
        print(f"  {k / kstar - 1:+10.1e} {r1:+11.3e} {r2:+11.3e} {r3:+11.3e} {runs['IC1']['pi'][i, j0]:+11.3e} "
              f"{runs['IC2']['pi'][i, j0]:+11.3e} {runs['IC1']['pi'][i, -1]:+12.3e} {runs['IC2']['pi'][i, -1]:+12.3e} "
              f"{den:+13.3e}")
    out.write_text(json.dumps({"alpha_M": m, "turn_on": ton, "C": C, "D1": D1, "D2": D2, "kstar": kstar,
                               "rows": rows}, indent=1, default=float))


if __name__ == "__main__":
    main()

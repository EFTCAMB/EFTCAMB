"""Verify fr-rgr-threshold-high-k-pk: time evolution at k = 1, 3, 10 /Mpc of designer f(R) (B0=1e-4) for the default
GR threshold (turn-on a=0.093) and threshold 1e-12 (turn-on a=0.01), via EFTCAMB.get_scale_evolution (harness evol_k).

Prints delta_c ratio (default / early-turn-on) and the effective mu (evol_mu) vs the QS mu (evol_QS_mu) around and
after the late turn-on, to see whether the late turn-on loses growth (ratio < 1) or injects a transient (ratio > 1).

Usage: python frthr_evol.py [run|show]
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

INI = harness.PARAMS_DIR / "4_Des_F_R_LCDM_2.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "frthr" / "evol"
KS = [1.0, 3.0, 10.0]
CASES = {"def": (INI, {}), "thr12": (INI, {"EFTCAMB_GR_threshold": "1e-12"}), "GR": (GR, {})}


def run():
    def job(item):
        name, (ini, ov) = item
        st = harness.run_model("A0", ini, OUT / f"{name}.npz", overrides=ov, evol_k=KS, threads=1, timeout=1500)
        return name, st["status"], st.get("elapsed")

    with ThreadPoolExecutor(3) as ex:
        for r in ex.map(job, CASES.items()):
            print(r)


def show():
    d = {n: harness.load(OUT / f"{n}.npz") for n in ("def", "thr12")}
    a = d["def"]["evol_a"]
    if a.ndim > 1:
        a = a[0]
    marks = [0.05, 0.09, 0.0931, 0.095, 0.1, 0.12, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0]
    idx = [int(np.argmin(abs(a - m))) for m in marks]
    for ik, k in enumerate(KS):
        print(f"k = {k} /Mpc")
        print(f"{'a':>8s} {'dc_def/dc_12':>13s} {'mu_def':>10s} "
              f"{'mu_12':>10s} {'QSmu':>10s} {'pi_def':>11s} {'pi_12':>11s}")
        for i in idx:
            cd, c12 = (d[n]["evol_clxc"][ik, i] for n in ("def", "thr12"))
            print(f"{a[i]:8.4f} {cd / c12:13.5f} {d['def']['evol_mu'][ik, i]:10.4f} "
                  f"{d['thr12']['evol_mu'][ik, i]:10.4f} {d['thr12']['evol_QS_mu'][ik, i]:10.4f} "
                  f"{d['def']['evol_pi'][ik, i]:11.3e} {d['thr12']['evol_pi'][ik, i]:11.3e}")


if __name__ == "__main__":
    run() if (len(sys.argv) > 1 and sys.argv[1] == "run") else show()

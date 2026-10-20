"""Diagnostic for metric_h_ic_violation.py: why is |h'_used - h'_00| O(1) in a ~ 1e-7..1e-5 even with the IC fix?

Runs (instrumented tree mh_instr, 5_hdsk_TG_1, turn-on 1e-8) evolve_metric_h = T with EFT_MH_ICFIX=1 and the default
F run, both logged, and prints the k ~ 10 and k ~ 1 /Mpc time series (accepted and trial derivs calls alike):
a, h'_used, h'_00, h''_ii, flags TC/RSA/no_phot, dphi.  For F also d h'_00/dtau (finite difference) vs h''_ii.

Usage: python metric_h_ic_diag.py
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

TREE = harness.AUDIT_ROOT / "work" / "verify" / "mh_instr"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "metric_h"
INI = harness.PARAMS_DIR / "5_hdsk_TG_1.ini"
COLS = "k tau a hdot hc hdd etak dgrho dphi dphidot rsa tc nophot dgq etakdot sigma".split()
CASES = {"Ticfix": ({"EFTCAMB_evolve_metric_h": "T"}, {"EFT_MH_ICFIX": "1"}),
         "F": ({"EFTCAMB_evolve_metric_h": "F"}, {})}


def run(item):
    name, (ov, env) = item
    binf = OUT / f"diag_{name}.bin"
    if not binf.exists():
        st = harness.run_model(TREE, INI, OUT / f"diag_{name}.npz", overrides=ov,
                               extra_env={**env, "EFT_MH_LOG": str(binf)}, threads=1, timeout=900)
        if st["status"] != "ok":
            raise RuntimeError(f"{name}: {st['status']}")
    raw = np.fromfile(binf, dtype=np.float64)
    raw = raw[: raw.size // 16 * 16].reshape(-1, 16)
    return name, {c: raw[:, i] for i, c in enumerate(COLS)}


with ThreadPoolExecutor(2) as ex:
    res = dict(ex.map(run, CASES.items()))

for kt in (10.0, 1.0):
    for name, d in res.items():
        ks = np.unique(d["k"])
        kk = ks[np.argmin(abs(np.log(ks / kt)))]
        m = d["k"] == kk
        o = np.argsort(d["tau"][m], kind="stable")
        x = {c: v[m][o] for c, v in d.items()}
        print(f"=== {name} k={kk:.4g}: {m.sum()} records; a range {x['a'].min():.2e}..{x['a'].max():.2e}")
        print(f"{'a':>10s} {'tau':>10s} {'hdot':>12s} {'hc':>12s} {'rel':>10s} {'hdd':>12s} {'TC':>3s} {'RSA':>3s} "
              f"{'nph':>3s} {'dphi':>10s} {'etak':>11s} {'dgrho':>11s}")
        for target in np.logspace(-8, -4, 17):
            i = int(np.argmin(abs(np.log(x["a"] / target))))
            print(f"{x['a'][i]:10.3e} {x['tau'][i]:10.3e} {x['hdot'][i]:+12.4e} {x['hc'][i]:+12.4e} "
                  f"{x['hdot'][i] / x['hc'][i] - 1:+10.2e} {x['hdd'][i]:+12.4e} {int(x['tc'][i]):3d} "
                  f"{int(x['rsa'][i]):3d} {int(x['nophot'][i]):3d} {x['dphi'][i]:+10.2e} {x['etak'][i]:+11.3e} "
                  f"{x['dgrho'][i]:+11.3e}")

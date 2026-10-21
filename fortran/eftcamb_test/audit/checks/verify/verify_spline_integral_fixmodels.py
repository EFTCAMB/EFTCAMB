"""Verifier (spline-integral-convention): rerun the spline-w model cases of verify_spline_integral_models.py on the
integral-only patched tree and compare with the EFTwDE = 1 references run on A0. If integral() is the whole root
cause, every case must now run and agree with the constant-w reference to quadrature/accuracy noise.
usage: /Users/marco/Software/miniforge3/bin/python verify_spline_integral_fixmodels.py
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
from verify_spline_integral_models import CASES, OUT  # noqa: E402

TREE = harness.AUDIT_ROOT / "work" / "verify" / "spline_integral" / "tree_intfix"
TESTS = {"DFR_ref": ["DFR_s9", "DFR_s9b", "DFR_s9def", "DFR_s10", "DFR_s10b"], "SS_ref": ["SS_s9b", "SS_s10b"]}


def run(name):
    ini, ov = CASES[name]
    out = OUT / f"fix_{name}.npz"
    js = out.with_suffix(".json")
    if js.exists():
        return name, json.loads(js.read_text())
    return name, harness.run_model(TREE, ini, out, overrides=ov, threads=1, timeout=900)


def main():
    names = [t for v in TESTS.values() for t in v] + ["DFR_ref", "SS_ref"]
    with ThreadPoolExecutor(4) as ex:
        sts = dict(ex.map(run, names))
    for ref, tests in TESTS.items():
        r = harness.load(OUT / f"{ref}.npz")  # A0 reference
        a = r["eft_a_grid"]
        c0 = harness.compare(OUT / f"fix_{ref}.npz", r)
        print(f"fix_{ref} (patched tree, EFTwDE=1) vs A0 {ref}: bit_identical={c0['bit_identical']} "
              f"worst Cl {c0['worst_cl']:.2e}")
        for t in tests:
            st = sts[t]
            if st["status"] != "ok":
                print(f"fix_{t:9s} {st['status']}  {st.get('log_tail', [''])[-1][:160]}")
                continue
            d = harness.load(OUT / f"fix_{t}.npz")
            c = harness.compare(d, r)
            g, g0 = d["eft_grhov_t"], r["eft_grhov_t"]
            m = (a > 1e-6) & (np.abs(g0) > 0)
            rel = np.max(np.abs(g[m] - g0[m]) / np.abs(g0[m]))
            print(f"fix_{t:9s} ok  max rel grhov_t diff {rel:.2e}  worst Cl {c['worst_cl']:.2e}  "
                  f"P(k) {c['worst_pk']:.2e}  sigma8 {c['sigma8']:.2e}  EFT fields worst "
                  f"{max(c['eft'].values()):.2e} ({max(c['eft'], key=c['eft'].get)})")


if __name__ == "__main__":
    main()

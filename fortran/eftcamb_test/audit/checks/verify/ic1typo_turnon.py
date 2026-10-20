"""Verifier (ic1-pidot-eftpid-typo): turn-on-time sensitivity test, independent of using the fixed build as truth.

The EFT_IC_type = 1 initial condition puts pi and pi' on the quasi-static attractor
    pi_att = -h0 E / (C + k^2 D1 + k^4 D2)
at the moment EFT perturbations are switched on (EFTCAMB_turn_on_time). If pi' is set consistently, the late-time
observables must not depend on the (early) turn-on time; a wrong pi' excites a homogeneous mode whose amplitude depends
on where it is injected. We run 5_Kmimic_1/2 with several turn-on times on A0 (typo) and FIX (D2 restored) and report
P(k)/P_GR - 1 at z = 0 and the spread across turn-on times.

Also checks that FIX is otherwise identical to A0: 1_EFT_GR on FIX vs A0 must be bit identical.
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

FIX = harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_pidfix"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1typo" / "turnon"
KS = (0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0)
TURNON = ("1.e-8", "3.e-9", "3.e-8", "1.e-7")
MODELS = sys.argv[1:] or ["5_Kmimic_1", "5_Kmimic_2"]


def job(args):
    tag, tree, model, over = args
    out = OUT / f"{tag}.npz"
    if not out.exists():
        st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=over, threads=1, timeout=900)
        if st["status"] != "ok":
            return tag, {"status": st["status"], "tail": st.get("log_tail", [])[-4:]}
    return tag, harness.load(out)


def pk_at(d, k):
    return float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))


def main():
    jobs = [("GR_A0", "A0", "1_EFT_GR", {}), ("GR_FIX", str(FIX), "1_EFT_GR", {})]
    for m in MODELS:
        for t in TURNON:
            for tn, tr in (("A0", "A0"), ("FIX", str(FIX))):
                jobs.append((f"{tn}_{m}_on{t}", tr, m, {"EFTCAMB_turn_on_time": t}))
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(job, jobs))
    summary = {}
    gr = res["GR_A0"]
    g2 = res["GR_FIX"]
    bit = all(np.array_equal(gr[k], g2[k]) for k in ("pk", "cl_total", "cl_lens_potential"))
    print(f"1_EFT_GR A0 vs FIX bit identical (pk, cl_total, cl_lens_potential): {bit}")
    summary["GR_bit_identical"] = bit
    for m in MODELS:
        print(f"\n{m}: P(k)/P_GR - 1 at z=0, k = " + ", ".join(f"{k:g}" for k in KS) + " /Mpc")
        for tn in ("A0", "FIX"):
            rows = []
            for t in TURNON:
                d = res[f"{tn}_{m}_on{t}"]
                if "status" in d:
                    print(f"  {tn:3s} turn_on={t:6s}  FAILED {d['status']} {d['tail']}")
                    continue
                r = [pk_at(d, k) / pk_at(gr, k) - 1 for k in KS]
                rows.append(r)
                print(f"  {tn:3s} turn_on={t:6s}  " + "  ".join(f"{x:+.4e}" for x in r))
            rows = np.array(rows)
            # spread of P(k) itself across turn-on times, relative: max/min - 1 of (1 + r)
            spread = (1 + rows).max(axis=0) / (1 + rows).min(axis=0) - 1
            print(f"  {tn:3s} spread max/min-1   " + "  ".join(f"{x:+.4e}" for x in spread))
            summary[f"{m}_{tn}"] = {"rows": rows.tolist(), "spread": spread.tolist()}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    main()

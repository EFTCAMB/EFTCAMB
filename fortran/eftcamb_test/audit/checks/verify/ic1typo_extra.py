"""Verifier (ic1-pidot-eftpid-typo): extra runs.

1. Other K-mouflage family inis (5_Kmimic_3, 5_Kmouflage_1..3) and a D2 != 0 model (5_Horava_eta_1): A0 vs FIX
   (FIX = A0 with only 09_EFTCAMB_IC.f90:589 corrected) -> worst Cl / P(k) change.
2. Convergence claim: 5_Kmimic_2 at accuracy_boost = 2 on A0 and FIX vs default; if the typo were the reason for
   non-convergence, FIX should converge while A0 does not.
3. C_L^phiphi and lensed TT change A0 vs FIX for 5_Kmimic_1/2 at default (from cached turn-on runs at 1e-8).
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
OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1typo" / "extra"
TON = harness.AUDIT_ROOT / "work" / "verify" / "ic1typo" / "turnon"
OTHER = ["5_Kmimic_3", "5_Kmouflage_1", "5_Kmouflage_2", "5_Kmouflage_3", "5_Horava_eta_1"]


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


def clmax(a, b, col, lo=2, hi=2500):
    x, y = a["cl_total"][lo:hi + 1, col], b["cl_total"][lo:hi + 1, col]
    r = np.abs(x / y - 1)
    return float(r.max()), int(np.argmax(r)) + lo


def main():
    jobs = []
    for m in OTHER:
        jobs += [(f"A0_{m}", "A0", m, {}), (f"FIX_{m}", str(FIX), m, {})]
    jobs += [("A0_5_Kmimic_2_boost2", "A0", "5_Kmimic_2", {"accuracy_boost": "2"}),
             ("FIX_5_Kmimic_2_boost2", str(FIX), "5_Kmimic_2", {"accuracy_boost": "2"})]
    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(job, jobs))
    summ = {}
    ks = (0.1, 1.0, 2.0, 5.0, 10.0, 20.0)
    print("1. A0 vs FIX, other models")
    for m in OTHER:
        a, f = res[f"A0_{m}"], res[f"FIX_{m}"]
        if "status" in a or "status" in f:
            print(f"  {m}: FAILED {a.get('status') if 'status' in a else ''} {f.get('status') if 'status' in f else ''}")
            continue
        c = harness.compare(a, f)
        pkr = [pk_at(a, k) / pk_at(f, k) - 1 for k in ks]
        print(f"  {m:16s} worst_cl={c['worst_cl']}  worst_pk={c['worst_pk']}  bit={c['bit_identical']}  "
              f"P_A0/P_FIX-1 @k={ks}: " + " ".join(f"{x:+.2e}" for x in pkr))
        summ[m] = {"worst_cl": c["worst_cl"], "worst_pk": c["worst_pk"], "pk_ratio": pkr}
    print("\n2. 5_Kmimic_2 convergence: default vs accuracy_boost=2 (lensed TT / EE max rel diff, ell 2..2500)")
    for tn in ("A0", "FIX"):
        d = harness.load(TON / f"{tn}_5_Kmimic_2_on1.e-8.npz")
        b = res[f"{tn}_5_Kmimic_2_boost2"]
        if "status" in b:
            print(f"  {tn}: boost2 FAILED {b['status']}")
            continue
        tt, ee = clmax(d, b, 0), clmax(d, b, 1)
        pkr = [pk_at(d, k) / pk_at(b, k) - 1 for k in ks]
        print(f"  {tn:3s} TT max {tt[0]:.3e} @ell {tt[1]}   EE max {ee[0]:.3e} @ell {ee[1]}   "
              f"P_def/P_boost2-1 @k={ks}: " + " ".join(f"{x:+.2e}" for x in pkr))
        summ[f"conv_{tn}"] = {"TT": tt, "EE": ee, "pk": pkr}
    print("\n3. lensing potential / lensed TT, A0 vs FIX at default (turn_on 1e-8)")
    for m in ("5_Kmimic_1", "5_Kmimic_2"):
        a = harness.load(TON / f"A0_{m}_on1.e-8.npz")
        f = harness.load(TON / f"FIX_{m}_on1.e-8.npz")
        pa, pf = a["cl_lens_potential"][:, 0], f["cl_lens_potential"][:, 0]
        ta, tf = a["cl_total"][:, 0], f["cl_total"][:, 0]
        ua, uf = a["cl_unlensed_scalar"][:, 0], f["cl_unlensed_scalar"][:, 0]
        print(f"  {m}: Cpp A0/FIX-1 @L=500,1000,2000,2500: "
              + " ".join(f"{pa[L] / pf[L] - 1:+.2e}" for L in (500, 1000, 2000, 2500))
              + "   lensed TT @ell=1000,2000,2500: " + " ".join(f"{ta[L] / tf[L] - 1:+.2e}" for L in (1000, 2000, 2500))
              + f"   unlensed TT max |rel| ell<=2500: {np.max(np.abs(ua[2:2501] / uf[2:2501] - 1)):.2e}")
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    main()

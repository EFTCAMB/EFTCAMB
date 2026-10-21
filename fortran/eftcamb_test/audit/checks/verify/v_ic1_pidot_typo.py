"""Verifier for ic1-pidot-eftpid-typo: independent IC-sensitivity check.
Runs 5_Kmimic_1 with (a) A0, (b) A0_icsw with pi'_ini = 0 (EFT_T5_PIDOT0=1, D function irrelevant),
(c) A0_icsw with centred-difference derivatives (EFT_T5_PIDER=1, fixed D), (d) A0 with turn-on time 1e-7.
If the fixed result is the attractor, (b),(c) agree with the fixed P(k) while A0 does not.
Also prints D1, D2 (from the EFT cache) to show k^4 D1 vs k^4 D2 magnitude."""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "t5-robust"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1_pidot_typo"
KS = (1.0, 2.0, 5.0, 10.0, 20.0)
M = "5_Kmimic_1"


def job(a):
    name, tree, model, kw = a
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, threads=1, timeout=600, **kw)
    if st["status"] != "ok":
        raise RuntimeError(f"{name}: {st['status']} {st.get('log_tail', [])[-3:]}")
    return name, harness.load(out)


def pk_at(d, k):
    return float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))


jobs = [("GR", "A0", "1_EFT_GR", {}), ("A0", "A0", M, {}), ("FIX", str(W / "A0_pidfix"), M, {}),
        ("PIDOT0", str(W / "A0_icsw"), M, {"extra_env": {"EFT_T5_PIDOT0": "1"}}),
        ("FIX_CDIFF", str(W / "A0_icsw"), M, {"extra_env": {"EFT_T5_PIDER": "1"}}),
        ("A0_ton1e-7", "A0", M, {"overrides": {"EFTCAMB_turn_on_time": "1.e-7"}}),
        ("FIX_ton1e-7", str(W / "A0_pidfix"), M, {"overrides": {"EFTCAMB_turn_on_time": "1.e-7"}})]
OUT.mkdir(parents=True, exist_ok=True)
with ThreadPoolExecutor(4) as ex:
    res = dict(ex.map(job, jobs))
gr = res["GR"]
print("P/P_GR-1 at k =", KS)
for n in [j[0] for j in jobs[1:]]:
    print(f"  {n:12s}", "  ".join(f"{pk_at(res[n], k) / pk_at(gr, k) - 1:+.4e}" for k in KS))
d = res["A0"]
a = d["eft_a_grid"]
for key in ("eft_EFTpiC", "eft_EFTpiD1", "eft_EFTpiD2", "eft_EFTpiE"):
    if key in d:
        v = d[key]
        i = np.argmin(abs(np.log(a) - np.log(1e-8)))
        print(key, "at a=1e-8..1e-6:", " ".join(f"{v[j]:+.3e}" for j in (i, i + 50, i + 100)))

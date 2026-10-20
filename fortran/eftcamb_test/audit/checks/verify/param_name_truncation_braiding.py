"""Verifier: RPH braiding as steplog function of Omega_DE. Full key RPHbraiding_ODE_delta vs truncated key
RPHbraiding_ODE_delt, two values each. Truncation predicts: full-key runs ignore the value; truncated-key runs respond."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness

OUT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
           "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/param-name-truncation")
OUT.mkdir(parents=True, exist_ok=True)
ini = harness.PARAMS_DIR / "3_RPH_KT_OmegaDE_1.ini"
base = {"RPHbraidingmodel_ODE": "8", "RPHbraiding_ODE_v1": "0.0", "RPHbraiding_ODE_v2": "0.05",
        "RPHbraiding_ODE_at": "0.5"}
runs = {"full_d0p3": {**base, "RPHbraiding_ODE_delta": "0.3"},
        "full_d0p6": {**base, "RPHbraiding_ODE_delta": "0.6"},
        "trunc_d0p3": {**base, "RPHbraiding_ODE_delt": "0.3"},
        "trunc_d0p6": {**base, "RPHbraiding_ODE_delt": "0.6"}}
res = {}
for name, ov in runs.items():
    st = harness.run_model("A0", ini, OUT / f"br_{name}.npz", overrides=ov, threads=1)
    tail = [t.strip() for t in (st.get("log_tail") or []) if "NaN" in t or "nstable" in t or "rror" in t][:2]
    print(name, st["status"], " / ".join(tail) if st["status"] != "ok" else "")
    res[name] = harness.load(OUT / f"br_{name}.npz") if st["status"] == "ok" else None
for a, b in [("full_d0p3", "full_d0p6"), ("trunc_d0p3", "trunc_d0p6")]:
    if res[a] is not None and res[b] is not None:
        c = harness.compare(res[a], res[b])
        print(f"{a} vs {b}: bit_identical={c['bit_identical']} worst_cl={c['worst_cl']:.3e}")

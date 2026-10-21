"""Verifier (param-name-truncation): model-level attribution with a fix tree (EFT_names_max_length 20 -> 64 only).

1. 3_RPH_KT_OmegaDE_1 + steplog kineticity(Omega_DE) (RPHkineticitymodel_ODE = 8):
   A0 documented key RPHkineticity_ODE_delta, A0 truncated key RPHkineticity_ODE_de, FIX documented key.
   Expected if the finding is right: A0-documented fails/differs, FIX-documented == A0-truncated.
2. 007p8 3+1 designer (audit ini 007p8-3p1.ini, sigma1_ODE linear with Designer3p1sigma1_ODE0 = 0.1):
   A0 vs FIX, and A0 with Designer3p1sigma1_ODE0 = 0 (should equal A0 with 0.1 if the key is never read).

usage: python param_trunc_model.py
"""

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "param_trunc"
FIX = OUT / "FIX"
INIS = HERE.parent.parent / "inis"

base = {"RPHkineticitymodel_ODE": "8", "RPHkineticity_ODE_v1": "0.01", "RPHkineticity_ODE_v2": "0.2",
        "RPHkineticity_ODE_at": "0.3"}
rph = harness.PARAMS_DIR / "3_RPH_KT_OmegaDE_1.ini"
p31 = INIS / "007p8-3p1.ini"
RUNS = {
    "rph_A0_doc": ("A0", rph, {**base, "RPHkineticity_ODE_delta": "0.5"}),
    "rph_A0_trunc": ("A0", rph, {**base, "RPHkineticity_ODE_de": "0.5"}),
    "rph_FIX_doc": (FIX, rph, {**base, "RPHkineticity_ODE_delta": "0.5"}),
    "rph_FIX_nodelta": (FIX, rph, dict(base)),
    "p31_A0": ("A0", p31, {}),
    "p31_A0_ode0": ("A0", p31, {"Designer3p1sigma1_ODE0": "0.0"}),
    "p31_FIX": (FIX, p31, {}),
    "p31_FIX_ode0": (FIX, p31, {"Designer3p1sigma1_ODE0": "0.0"}),
}


def main():
    st = {}
    for name, (tree, ini, ov) in RUNS.items():
        s = harness.run_model(tree, ini, OUT / f"{name}.npz", overrides=ov, threads=1, timeout=900)
        tail = [t.strip() for t in (s.get("log_tail") or []) if "NaN" in t or "unstable" in t or "failed" in t]
        st[name] = s["status"]
        print(f"{name:18s} {s['status']:12s} {' / '.join(tail[:2]) if s['status'] != 'ok' else ''}", flush=True)

    def cmp(a, b):
        if st.get(a) != "ok" or st.get(b) != "ok":
            print(f"compare {a} vs {b}: n/a ({st.get(a)}, {st.get(b)})")
            return
        c = harness.compare(OUT / f"{a}.npz", OUT / f"{b}.npz")
        eft = max(c.get("eft", {}).values() or [0.0])
        print(f"compare {a} vs {b}: worst_cl {c['worst_cl']:.3e} worst_pk {c['worst_pk']:.3e} "
              f"sigma8 {c.get('sigma8', float('nan')):.3e} eft_max {eft:.3e} bit_identical {c['bit_identical']}")

    cmp("rph_FIX_doc", "rph_A0_trunc")
    cmp("rph_FIX_doc", "rph_FIX_nodelta")
    cmp("p31_A0", "p31_A0_ode0")
    cmp("p31_FIX", "p31_FIX_ode0")
    cmp("p31_FIX_ode0", "p31_A0")
    for n in ("rph_A0_trunc", "rph_FIX_doc", "p31_A0", "p31_FIX"):
        f = OUT / f"{n}.npz"
        if f.exists():
            m = harness.load(f)["meta"]
            ep = m.get("eft_params") or {}
            sel = {k: v for k, v in ep.items() if "delta" in k or "sigma1_ODE" in k or "_de" in k}
            print(f"{n}: eft_params subset {json.dumps(sel)}")


if __name__ == "__main__":
    main()

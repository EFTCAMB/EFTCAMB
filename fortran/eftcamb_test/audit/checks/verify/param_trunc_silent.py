"""Verifier (param-name-truncation): a silent wrong-physics case (no error raised).

3_RPH_KT_OmegaDE_1 with the kineticity given as a Fourier series in Omega_DE (RPHkineticitymodel_ODE = 7):
alpha_K = a0 + a1 cos(2 pi (Omega_DE - phase)), a0 = 0.1, a1 = 0.05, phase = RPHkineticity_ODE_phase (23 chars)
(alpha_K >= 0.05 > 0 for any phase, so the ghost condition holds).
If the key is truncated, setting the phase has no effect on A0 (bit-identical to no phase), while the fix tree
(EFT_names_max_length = 64 only) responds to it; without the phase key A0 and FIX must be bit-identical.

usage: python param_trunc_silent.py
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "param_trunc"
FIX = OUT / "FIX"
ini = harness.PARAMS_DIR / "3_RPH_KT_OmegaDE_1.ini"
base = {"RPHkineticitymodel_ODE": "7", "RPHkineticity_ODE_Fourier_order_cos": "1", "RPHkineticity_ODE_a0": "0.1",
        "RPHkineticity_ODE_a1": "0.05"}
PH = {"RPHkineticity_ODE_phase": "0.25"}
RUNS = {
    "kt_A0_phase": ("A0", {**base, **PH}),
    "kt_A0_nophase": ("A0", base),
    "kt_FIX_phase": (FIX, {**base, **PH}),
    "kt_FIX_nophase": (FIX, base),
    "kt_A0_truncphase": ("A0", {**base, "RPHkineticity_ODE_ph": "0.25"}),
}


def main():
    st = {}
    for name, (tree, ov) in RUNS.items():
        s = harness.run_model(tree, ini, OUT / f"{name}.npz", overrides=ov, threads=1, timeout=900)
        st[name] = s["status"]
        tail = [t.strip() for t in (s.get("log_tail") or []) if "NaN" in t or "unstable" in t or "failed" in t]
        print(f"{name:18s} {s['status']:12s} {' / '.join(tail[:2]) if s['status'] != 'ok' else ''}", flush=True)

    def cmp(a, b):
        if st.get(a) != "ok" or st.get(b) != "ok":
            print(f"compare {a} vs {b}: n/a ({st.get(a)}, {st.get(b)})")
            return
        c = harness.compare(OUT / f"{a}.npz", OUT / f"{b}.npz")
        eft = c.get("eft", {})
        ab = eft.get("EFTalphaB", eft.get("alphaB", max(eft.values() or [0.0])))
        print(f"compare {a:16s} vs {b:16s}: worst_cl {c['worst_cl']:.3e} worst_pk {c['worst_pk']:.3e} "
              f"sigma8 {c.get('sigma8', float('nan')):.3e} eft_max {max(eft.values() or [0.0]):.3e} "
              f"bit_identical {c['bit_identical']}")

    cmp("kt_A0_phase", "kt_A0_nophase")
    cmp("kt_FIX_nophase", "kt_A0_nophase")
    cmp("kt_FIX_phase", "kt_FIX_nophase")
    cmp("kt_FIX_phase", "kt_A0_truncphase")


if __name__ == "__main__":
    main()

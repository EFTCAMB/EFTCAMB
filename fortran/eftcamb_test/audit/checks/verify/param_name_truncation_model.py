"""Verifier: full 007p8 3+1 designer runs. If the 21-char key Designer3p1sigma1_ODE0 is truncated, setting it to 0.1
must give the same result as 0.0, while the truncated key Designer3p1sigma1_OD = 0.1 changes the result."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness

OUT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
           "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/param-name-truncation")
OUT.mkdir(parents=True, exist_ok=True)
ini = Path(__file__).resolve().parents[2] / "inis" / "verify_007p8_trunc.ini"
runs = {"full_0p1": {"Designer3p1sigma1_ODE0": "0.1"},
        "full_0p0": {"Designer3p1sigma1_ODE0": "0.0"},
        "trunc_0p1": {"Designer3p1sigma1_ODE0": "0.0", "Designer3p1sigma1_OD": "0.1"},
        "full_0p3": {"Designer3p1sigma1_ODE0": "0.3"}}
res = {}
for name, ov in runs.items():
    st = harness.run_model("A0", ini, OUT / f"{name}.npz", overrides=ov, threads=1)
    print(name, st["status"])
    res[name] = harness.load(OUT / f"{name}.npz") if st["status"] == "ok" else None
for a, b in [("full_0p1", "full_0p0"), ("full_0p1", "full_0p3"), ("trunc_0p1", "full_0p0")]:
    if res[a] is not None and res[b] is not None:
        c = harness.compare(res[a], res[b])
        print(f"{a} vs {b}: bit_identical={c['bit_identical']} worst_cl={c['worst_cl']:.3e} "
              f"eft_max={max(c['eft'].values()) if c['eft'] else 0:.3e}")

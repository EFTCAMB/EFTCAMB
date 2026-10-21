"""Verifier (t3-pi-ic-singular-denominator): attribute P(k) spikes of 3_RPH_MK_const_1 to the EFT_IC_type=1 pi IC.
Runs P(k) on the knob tree (see t3_piic_build_knob.py) with AUDIT_VPI_REG unset / 0.1 / 1e30 (pi=pi'=0 always),
and a dense-k scan of pi(a_ini), delta_cdm(a=1) across k* = sqrt(-C/D1) on A0."""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

KNOB = harness.AUDIT_ROOT / "work" / "verify" / "t3piic_knob"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "t3piic"
INI = harness.PARAMS_DIR / "3_RPH_MK_const_1.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
PK = {"get_scalar_cls": "F", "get_tensor_cls": "F", "do_lensing": "F", "transfer_kmax": "10",
      "num_redshiftwindows": "0"}


def job(a):
    name, tree, ini, ov, env, ps = a
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, ini, out, overrides={**PK, **ov}, python_set=ps, extra_env=env, threads=1,
                           timeout=600)
    return name, st, out


jobs = [("GR", "A0", GR, {}, {}, {"Accuracy.IntTolBoost": 2.0})]
for m in ["1e-07", "1e-06"]:
    jobs.append((f"A0_m{m}", "A0", INI, {"RPHmassP0": m}, {}, {}))
    for lab, f in [("unset", None), ("reg0.1", "0.1"), ("zero", "1e30")]:
        jobs.append((f"K_{lab}_m{m}", KNOB, INI, {"RPHmassP0": m}, {} if f is None else {"AUDIT_VPI_REG": f}, {}))
    # shift turn-on time: k* should move
    jobs.append((f"A0_ton0.012_m{m}", "A0", INI, {"RPHmassP0": m, "EFTCAMB_turn_on_time": "0.012"}, {}, {}))
with ThreadPoolExecutor(4) as ex:
    res = {n: (st, o) for n, st, o in ex.map(job, jobs)}
ref = harness.load(res["GR"][1])
k = ref["pk_k"]
for n, (st, o) in res.items():
    if n == "GR":
        continue
    if st["status"] != "ok":
        print(n, st["status"]); continue
    d = harness.load(o)
    r = np.abs(d["pk"][0] / ref["pk"][0] - 1)
    i = int(np.argmax(r))
    print(f"{n:22s} max|P/P_GR-1|={r[i]:.2e} at k={k[i]:.4g}  median={np.median(r):.1e}")
a0 = harness.load(res["A0_m1e-06"][1])["pk"][0]
for lab in ["unset", "reg0.1", "zero"]:
    b = harness.load(res[f"K_{lab}_m1e-06"][1])["pk"][0]
    print(f"knob {lab} vs A0 (m=1e-6): max rel diff {np.max(np.abs(b / a0 - 1)):.2e}")

# dense k scan across k* on A0 (alpha_M=1e-6)
ks = [float(x) for x in np.linspace(0.160, 0.180, 21)]
o = OUT / "evol_dense.npz"
st = harness.run_model("A0", INI, o, overrides={"RPHmassP0": "1e-06", "get_tensor_cls": "F", "do_lensing": "F",
                       "l_max_scalar": "100", "num_redshiftwindows": "0"}, evol_k=ks, threads=1, timeout=600)
print("evol", st["status"])
d = harness.load(o)
a = np.atleast_2d(d["evol_a"])[0]
j = int(np.searchsorted(a, 0.01))
print(f"a_ini sample = {a[j]:.5f}; fields: {[x for x in d if x.startswith('evol_')][:60]}")
for n, kv in enumerate(ks):
    C, D1, D2, E = (d[f"evol_{f}"][n, j] for f in ("EFTpiC", "EFTpiD1", "EFTpiD2", "EFTpiE"))
    den = C + kv**2 * D1 + kv**4 * D2
    pi0 = d["evol_pi"][n, j]
    pi1 = d["evol_pi"][n, -1]
    dc = d["evol_clxc"][n, -1] if "evol_clxc" in d else np.nan
    print(f"k={kv:.4f} den={den:+.2e} -h0E/den(*)={-E/den:+.2e} pi(a_ini)={pi0:+.3e} pi(1)={pi1:+.3e} clxc(1)={dc:+.6e}")

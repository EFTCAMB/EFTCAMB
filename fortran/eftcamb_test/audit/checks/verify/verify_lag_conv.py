"""Verifier check for fullmap-bg-eft-stale-background, independent of t2d-static's LAG tree.

Tree work/verify/lag/CONV = A0 + patch_conv_lag.py (derivs iterates the full-model background to its fixed point on
a copy of the cache and logs CONVMAX = max over derivs calls (a >= 1e-3) of the used-vs-fixed-point mismatch,
normalised by natural scales; EFT_CONV_FIX=1 continues with the fixed-point cache).
Runs (Beyond Horndeski, audit ini 008p5-beyond-horndeski.ini, EFTCAMB_skip_stability = T, 1 thread):
  a0      A0, default accuracy           conv     CONV, used values kept      convfix   CONV, fixed point used
  a0_acc2 A0, accuracy_boost=2           conv_acc2                            convfix_acc2
Compares spectra: lag effect (convfix vs a0) at accuracy 1 and 2, accuracy noise (a0_acc2 vs a0), control (conv vs a0).

Usage: python verify_lag_conv.py [model_ini]   (outputs in work/verify/lag/runs/)
"""

from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
CONV = AUDIT / "work" / "verify" / "lag" / "CONV"
OUT = AUDIT / "work" / "verify" / "lag" / "runs"
INI = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parents[1] / "inis" / "008p5-beyond-horndeski.ini"
NAMES = {1: "H", 2: "Hdot", 3: "EFTc", 4: "EFTLambda", 5: "EFTcdot", 6: "EFTLambdadot", 7: "EFTc(plain rel)",
         8: "Hdotdot"}
SCALE = {1: "|H_u/H_fp-1|", 2: "/H^2", 3: "/(3H^2)", 4: "/(3H^2)", 5: "/(3H^3)", 6: "/(3H^3)", 7: "/|EFTc_fp|",
         8: "/H^3"}


def jobs(tag: str):
    base = {"EFTCAMB_skip_stability": "T"}
    acc2 = dict(base, accuracy_boost="2")
    return [
        (f"{tag}_a0", "A0", base, {}),
        (f"{tag}_conv", CONV, base, {}),
        (f"{tag}_convfix", CONV, base, {"EFT_CONV_FIX": "1"}),
        (f"{tag}_a0_acc2", "A0", acc2, {}),
        (f"{tag}_conv_acc2", CONV, acc2, {}),
        (f"{tag}_convfix_acc2", CONV, acc2, {"EFT_CONV_FIX": "1"}),
    ]


def convmax(log: Path) -> dict:
    best: dict = {}
    summ = None
    for line in log.read_text(errors="replace").splitlines():
        if line.startswith("CONVMAX"):
            v = line.split()
            i, it = int(v[1]), int(v[2])
            a, k, tau, used, fp, m, H = map(float, v[3:10])
            best[NAMES[i]] = dict(metric=m, scale=SCALE[i], a=a, k=k, tau=tau, used=used, fixed_point=fp, iters=it)
        elif line.startswith("CONVSUM"):
            summ = line
    return {"max": best, "last_sum": summ}


def worst(cmp: dict) -> dict:
    cl = cmp["cl"]
    r = {}
    for s in ("unlensed_scalar", "lens_potential", "total"):
        if s in cl:
            r[s] = {k: float(v) for k, v in cl[s].items()}
    r["pk_z0"] = float(cmp["pk"][0]) if cmp.get("pk") else None
    r["sigma8"] = cmp.get("sigma8")
    return r


def pk_detail(a: dict, b: dict) -> dict:
    """max rel diff of P(k, z=0) for k ranges."""
    if a["pk_k"].shape == b["pk_k"].shape and np.allclose(a["pk_k"], b["pk_k"], rtol=1e-12):
        k, pa, pb = a["pk_k"], a["pk"][0], b["pk"][0]
    else:
        k = np.logspace(np.log10(max(a["pk_k"][0], b["pk_k"][0])), np.log10(min(a["pk_k"][-1], b["pk_k"][-1])), 400)
        pa = np.exp(np.interp(np.log(k), np.log(a["pk_k"]), np.log(a["pk"][0])))
        pb = np.exp(np.interp(np.log(k), np.log(b["pk_k"]), np.log(b["pk"][0])))
    d = np.abs(pb / pa - 1)
    out = {}
    for lo, hi in [(0, 0.1), (0.1, 1), (1, 5), (5, 100)]:
        m = (k >= lo) & (k < hi)
        if m.any():
            j = np.argmax(d[m])
            out[f"k[{lo},{hi})"] = (float(d[m][j]), float(k[m][j]))
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    tag = INI.stem
    js = jobs(tag)

    def run(j):
        name, tree, ov, env = j
        npz = OUT / f"{name}.npz"
        if npz.exists() and json.loads((OUT / f"{name}.json").read_text()).get("status") == "ok":
            return name, {"status": "ok", "cached": True}
        return name, harness.run_model(tree, INI, npz, overrides=ov, threads=1, timeout=3000, extra_env=env)

    with ThreadPoolExecutor(4) as ex:
        st = dict(ex.map(run, js))
    res = {"ini": str(INI), "status": {k: v["status"] for k, v in st.items()}}
    print(json.dumps(res["status"], indent=1))
    ok = {k for k, v in st.items() if v["status"] == "ok"}
    for k in (f"{tag}_conv", f"{tag}_convfix", f"{tag}_conv_acc2", f"{tag}_convfix_acc2"):
        log = OUT / f"{k}.log"
        if log.exists():
            res[f"lag_{k}"] = convmax(log)
    pairs = {
        "control conv vs a0 (recompile only)": (f"{tag}_a0", f"{tag}_conv"),
        "LAG EFFECT convfix vs a0 (acc 1)": (f"{tag}_a0", f"{tag}_convfix"),
        "LAG EFFECT convfix_acc2 vs a0_acc2 (acc 2)": (f"{tag}_a0_acc2", f"{tag}_convfix_acc2"),
        "accuracy noise a0_acc2 vs a0": (f"{tag}_a0", f"{tag}_a0_acc2"),
        "accuracy noise convfix_acc2 vs convfix": (f"{tag}_convfix", f"{tag}_convfix_acc2"),
        "a0 (acc1) vs convfix_acc2 (best)": (f"{tag}_convfix_acc2", f"{tag}_a0"),
        "control conv_acc2 vs a0_acc2": (f"{tag}_a0_acc2", f"{tag}_conv_acc2"),
    }
    res["compare"] = {}
    for lab, (x, y) in pairs.items():
        if x in ok and y in ok:
            A, B = harness.load(OUT / f"{x}.npz"), harness.load(OUT / f"{y}.npz")
            c = harness.compare(A, B, lmax=2500)
            res["compare"][lab] = {"spectra": worst(c), "pk_bins": pk_detail(A, B)}
    for k, v in res.items():
        if k.startswith("lag_"):
            print(f"\n{k}: {v['last_sum']}")
            for f, d in v["max"].items():
                print(f"  {f:16s} {d['scale']:14s} max = {d['metric']:.3e} at a={d['a']:.4f} k={d['k']:.3e} "
                      f"used={d['used']:+.5e} fp={d['fixed_point']:+.5e} iters={d['iters']}")
    for lab, d in res["compare"].items():
        s = d["spectra"]
        us = s.get("unlensed_scalar", {})
        lp = s.get("lens_potential", {})
        print(f"\n{lab}:\n  unlensed TT {us.get('TT', float('nan')):.3e} EE {us.get('EE', float('nan')):.3e} "
              f"TE {us.get('TE', float('nan')):.3e}  phiphi {lp.get('PP', float('nan')):.3e}  P(k,z=0) {s['pk_z0']}  "
              f"sigma8 {s['sigma8']}")
        print("  P(k) by k bin (maxrel, k):", {k: (f"{v[0]:.2e}", f"{v[1]:.2f}") for k, v in d["pk_bins"].items()})
    (OUT / f"{tag}_summary.json").write_text(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()

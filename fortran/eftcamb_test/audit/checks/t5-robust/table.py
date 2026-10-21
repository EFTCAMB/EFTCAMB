"""Track 5: table of every variation run vs the A0 default run, with the noise floor of the same model.

noise(model) = max over the available high-accuracy reruns of the same model:
    AB2  = accuracy_boost 2                       (work/t5-robust/runs/AB2, this track)
    ITB2 = Accuracy.IntTolBoost 2                 (work/t10a-regression/runs_A0_itb2, track t10a)
    AB2I = accuracy_boost 2 + IntTolBoost 2       (work/t10a-regression/runs_A0_ab2, track t10a)
Observables: scal = worst lensed/unlensed scalar C_l (TT, EE, BB, TE/sqrt(TT EE), PP), ell 2..2500, max relative;
TT30 = lensed TT for ell >= 30; pk = P(k, z=0) max relative; s8 = sigma8 relative.
A variation is FLAGged when an observable changes by more than max(1e-4, noise of that observable).

Usage: python table.py VAR[,VAR...] [--models m1,m2] [--json out.json]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import analyse  # noqa: E402
import harness  # noqa: E402

R = harness.AUDIT_ROOT
WORK = R / "work" / "t5-robust"
NOISE_DIRS = {"AB2": WORK / "runs" / "AB2", "ITB2": R / "work" / "t10a-regression" / "runs_A0_itb2",
              "AB2I": R / "work" / "t10a-regression" / "runs_A0_ab2"}
OBS = ("scal", "TT30", "pk", "s8")


def ref_path(model: str) -> Path:
    p = R / "runs" / "A0" / f"{model}.npz"
    return p if p.exists() else R / "work" / "prep" / "runs" / f"{model}.npz"


def obs(A: dict, B: dict) -> dict:
    m = analyse.metrics(A, B)
    out = {"scal": m["worst_scalar_cl"], "pk": m.get("pk", float("nan")), "s8": m.get("sigma8", float("nan")),
           "TT30": analyse.frac(A["cl_total"][:, 0], B["cl_total"][:, 0], lmin=30),
           "ten": m["worst_tensor_cl"]}
    out["bitid"] = bool(np.array_equal(A["cl_total"], B["cl_total"]) and np.array_equal(A["pk"], B["pk"]))
    return out


_cache: dict = {}


def load(p: Path) -> dict:
    if p not in _cache:
        _cache[p] = harness.load(p)
    return _cache[p]


def noise(model: str) -> dict:
    out = {}
    for name, d in NOISE_DIRS.items():
        f = d / f"{model}.npz"
        js = f.with_suffix(".json")
        if f.exists() and (not js.exists() or json.loads(js.read_text())["status"] == "ok"):
            out[name] = obs(load(f), load(ref_path(model)))
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("variations")
    p.add_argument("--models", default="")
    p.add_argument("--json", default="")
    a = p.parse_args()
    allres = {}
    for var in a.variations.split(","):
        d = WORK / "runs" / var
        print(f"=== {var}")
        for js in sorted(d.glob("*.json")):
            model = js.stem
            if a.models and model not in a.models.split(","):
                continue
            st = json.loads(js.read_text())
            if st["status"] != "ok":
                print(f"  {model:30s} {st['status']}  " + " | ".join(x.strip() for x in st.get("log_tail", [])[-3:])[:200])
                allres.setdefault(var, {})[model] = {"status": st["status"]}
                continue
            v = obs(load(js.with_suffix(".npz")), load(ref_path(model)))
            nz = noise(model)
            thr = {o: max([1e-4] + [n[o] for n in nz.values() if np.isfinite(n[o])]) for o in OBS}
            flags = [o for o in OBS if not v[o] <= thr[o]]
            allres.setdefault(var, {})[model] = {"status": "ok", **v, "noise": nz, "thr": thr, "flags": flags}
            ns = " ".join(f"{k}:{n['scal']:.0e}/{n['pk']:.0e}" for k, n in nz.items())
            print(f"  {model:30s} {'BITID' if v['bitid'] else '     '} scal={v['scal']:.2e} TT30={v['TT30']:.2e} "
                  f"pk={v['pk']:.2e} s8={v['s8']:.2e} | noise(scal/pk) {ns:34s} {'FLAG ' + ','.join(flags) if flags else ''}")
    if a.json:
        Path(a.json).write_text(json.dumps(allres, indent=1, default=float))


if __name__ == "__main__":
    main()

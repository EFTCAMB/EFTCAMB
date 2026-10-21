"""Verify fr-rgr-threshold-high-k-pk: designer f(R) (4_Des_F_R_LCDM_2, B0=1e-4) P(k) vs EFTCAMB_GR_threshold.

Runs through harness.run_model (machine-wide slots).  Each config is NAME:key=val;key=val (overrides of the model ini).
Keys ending in "@env" are passed as environment variables instead (e.g. EFT_T5_PIDOT0@env=1).
A GR reference (1_EFT_GR.ini) with the same common overrides is run too, so P_fR/P_GR - 1 can be printed per config.

Usage: python frthr_scan.py TREE OUTSUBDIR "common_k=v;..." CONFIG [CONFIG ...]
       (common may be "" ; use jobs via env FRTHR_JOBS, default 4)
"""

import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

INI = harness.PARAMS_DIR / "4_Des_F_R_LCDM_2.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
KS = (0.1, 0.3, 1.0, 3.0, 10.0)


def parse(s: str) -> dict:
    return dict(kv.split("=", 1) for kv in s.split(";") if kv.strip())


def main() -> None:
    tree, sub, common = sys.argv[1], sys.argv[2], parse(sys.argv[3])
    configs = [(c.split(":", 1)[0], parse(c.split(":", 1)[1] if ":" in c else "")) for c in sys.argv[4:]]
    out = harness.AUDIT_ROOT / "work" / "verify" / "frthr" / sub
    jobs = int(os.environ.get("FRTHR_JOBS", "4"))

    def job(item):
        name, ov = item
        ini = GR if name == "GR" else INI
        npz = out / f"{name}.npz"
        if npz.exists() and npz.with_suffix(".json").exists():
            return name, {"status": "ok", "elapsed": 0.0, "cached": True}, npz
        allov = {**common, **ov}
        env = {k[:-4]: v for k, v in allov.items() if k.endswith("@env")}
        allov = {k: v for k, v in allov.items() if not k.endswith("@env")}
        st = harness.run_model(tree, ini, npz, overrides=allov, threads=1, timeout=1500, extra_env=env)
        return name, st, npz

    items = [("GR", {})] + configs
    with ThreadPoolExecutor(jobs) as ex:
        res = list(ex.map(job, items))
    ref = harness.load(out / "GR.npz")
    k = ref["pk_k"]
    ik = [int(np.argmin(abs(k - x))) for x in KS]
    base = None
    print(f"tree={tree} common={common}  kmax={k.max():.3g}")
    print(f"{'config':>14s} {'status':>7s} {'t[s]':>6s} {'a_on':>8s} {'sigma8':>9s}  P_fR/P_GR-1 at "
          + " ".join(f"k={k[i]:.3g}" for i in ik))
    rows = {}
    for name, st, npz in res:
        if st["status"] != "ok":
            print(f"{name:>14s} {st['status']:>7s}")
            continue
        d = harness.load(npz)
        log = npz.with_suffix(".log").read_text()
        rgr = re.findall(r"Return to GR time: *([-+0-9.E]+)", log)
        a_on = float(rgr[-1]) if rgr else float("nan")
        assert np.allclose(d["pk_k"], k), "k grids differ"
        r = d["pk"][0] / ref["pk"][0] - 1
        rows[name] = d["pk"][0]
        print(f"{name:>14s} {st['status']:>7s} {st.get('elapsed', 0):6.1f} {a_on:8.4f} {d['sigma8'][0]:9.6f}  "
              + " ".join(f"{r[i]:+10.4e}" for i in ik))
    base = configs[0][0]
    if base in rows:
        print(f"relative to {base}: P/P_{base}-1")
        for name in rows:
            if name in ("GR", base):
                continue
            r = rows[name] / rows[base] - 1
            print(f"{name:>14s} " + " ".join(f"{r[i]:+10.4e}" for i in ik))


if __name__ == "__main__":
    main()

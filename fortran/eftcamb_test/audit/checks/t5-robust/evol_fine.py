"""Track 5: pi(a), pi'(a) right after the EFT perturbation turn-on on a fine scale-factor grid (IC1 attractor check).

For each model: run camb.get_background on the chosen tree, read the effective turn-on (after the RGR shift) from
results.Params, then call EFTCAMB.get_scale_evolution on a = a_on * logspace(0, 1, NA) for several k and compare pi
with the quasi-static solution pi_QS = -h0 E / (C + k^2 D1 + k^4 D2) built from the same cache.

Diagnostics per k:
  r0       pi/pi_QS at the first grid point after the turn-on (IC1 sets r0 = 1 by construction)
  dev_early  max |pi/pi_QS - 1| over a in [a_on, 1.2 a_on]
  dev_late   max |pi/pi_QS - 1| over a in [3 a_on, 10 a_on]
  osc      number of sign changes of d(pi/pi_QS)/dln a over a in [a_on, 3 a_on] (transient oscillations)
Usage: python evol_fine.py MODEL[,MODEL...] [--tree A0] [--set k=v ...] [--k 0.01,0.1,1,10] [--json out]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t5-robust"
CHILD = r"""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
ks = [float(x) for x in sys.argv[3].split(',')]
na = int(sys.argv[4])
res = camb.get_background(cp)
a_on = float(res.Params.EFTCAMB.EFTCAMB_pert_turn_on)
a = a_on * np.logspace(0.0, 1.0, na)
a = a[a < 1.0]
a[0] *= 1.0 + 1e-6
fields, arr = cp.EFTCAMB.get_scale_evolution(res, ks, a)
out = {'a_on': a_on, 'a': arr['a'][0].tolist(), 'h0': cp.H0 / 299792.458}
for f in ('pi', 'pidot', 'EFTpiC', 'EFTpiD1', 'EFTpiD2', 'EFTpiE', 'EFTpiA1', 'EFTpiA2', 'EFTpiB1', 'EFTpiB2', 'k'):
    out[f] = arr[f].tolist()
print('JSON' + json.dumps(out))
"""


def run(tree, model, overrides, ks, na=400):
    ini = harness.PARAMS_DIR / f"{model}.ini"
    if not ini.exists():
        ini = Path(model)
    tmp = tempfile.mkdtemp(dir=WORK / "evolfine_tmp")
    w = os.path.join(tmp, "w.ini")
    lines = [f"DEFAULT({ini.resolve()})", f"output_root = {tmp}/r", "EFTCAMB_write_background = F",
             "feedback_level = 0"] + [f"{k} = {v}" for k, v in overrides.items()]
    Path(w).write_text("\n".join(lines) + "\n")
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", CHILD, str(harness.tree_path(tree)), w, ",".join(map(str, ks)), str(na)],
                       capture_output=True, text=True, cwd=tmp, env=env, timeout=3000)
    for line in p.stdout.splitlines():
        if line.startswith("JSON"):
            return json.loads(line[4:])
    return {"error": f"rc={p.returncode}: " + (p.stderr or p.stdout)[-500:]}


def diagnose(d, ks):
    a = np.array(d["a"])
    a_on = d["a_on"]
    out = {"a_on": a_on}
    for ik, k in enumerate(ks):
        pi = np.array(d["pi"][ik])
        C, D1, D2, E = (np.array(d[f][ik]) for f in ("EFTpiC", "EFTpiD1", "EFTpiD2", "EFTpiE"))
        den = C + k**2 * D1 + k**4 * D2
        with np.errstate(divide="ignore", invalid="ignore"):
            qs = np.where(den != 0, -d["h0"] * E / den, np.nan)
            r = pi / qs
        good = np.isfinite(r)
        if not good.any():
            out[str(k)] = {"note": "pi_QS undefined (C + k^2 D = 0)", "max_pi": float(np.max(np.abs(pi)))}
            continue
        early = good & (a <= 1.2 * a_on)
        late = good & (a >= 3 * a_on) & (a <= 10 * a_on)
        win = good & (a <= 3 * a_on)
        dr = np.diff(r[win])
        osc = int(np.sum(np.sign(dr[1:]) != np.sign(dr[:-1]))) if dr.size > 2 else 0
        out[str(k)] = {
            "r0": float(r[good][0]),
            "dev_early": float(np.nanmax(np.abs(r[early] - 1))) if early.any() else float("nan"),
            "dev_late": float(np.nanmax(np.abs(r[late] - 1))) if late.any() else float("nan"),
            "osc": osc,
        }
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("models")
    p.add_argument("--tree", default="A0")
    p.add_argument("--set", nargs="*", default=[])
    p.add_argument("--k", default="0.01,0.1,1,10")
    p.add_argument("--na", type=int, default=400)
    p.add_argument("--json", default="")
    a = p.parse_args()
    (WORK / "evolfine_tmp").mkdir(parents=True, exist_ok=True)
    ks = [float(x) for x in a.k.split(",")]
    ov = dict(s.split("=", 1) for s in a.set)
    allres = {}
    if a.json and Path(a.json).exists():
        allres = json.loads(Path(a.json).read_text())
    for m in a.models.split(","):
        if m in allres:
            continue
        d = run(a.tree, m, ov, ks, a.na)
        if "error" in d:
            res = {"error": d["error"][-300:]}
        else:
            res = diagnose(d, ks)
            raw = WORK / "evolfine" / f"{Path(m).stem}{'_' + a.tree if a.tree != 'A0' else ''}.json"
            raw.parent.mkdir(exist_ok=True)
            raw.write_text(json.dumps(d))
        allres[m] = res
        print(m, json.dumps(res), flush=True)
        if a.json:
            Path(a.json).write_text(json.dumps(allres, indent=1))


if __name__ == "__main__":
    main()

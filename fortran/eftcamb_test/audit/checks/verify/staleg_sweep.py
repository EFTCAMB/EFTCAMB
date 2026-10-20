"""Breadth sweep for ic23-stale-cache-gammas on the frozen A0 release build via get_time_evolution (1 thread).

For each model and EFT_IC_type in (2, 3): q = [k, k, k], output delta_photon just after taustart (adiabatic IC value),
last requested time a_last = 0.05 or 1.  A correct `initial` gives the same IC for every mode and every a_last.
Reports: mode-1 finiteness, |dg(mode2, a_last=0.05)/dg(mode2, a_last=1) - 1| (stale dependence on the previous mode),
and for IC 2 the Omega0 inverted from dg vs EFTOmegaV(a_eval) from the A0 reference run (get_eft_functions).
usage: python staleg_sweep.py model [model ...]   (or 'all')"""

import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

TREE = str(harness.TREES["A0"])
OUT = harness.AUDIT_ROOT / "work" / "verify" / "staleg" / "sweep"
OUT.mkdir(parents=True, exist_ok=True)
K = 0.05
CHILD = r"""
import json, sys
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_background(pars)
k = float(sys.argv[3])
tst = 0.001 / k
out = {"a_eval": max(pars.EFTCAMB.EFTCAMB_pert_turn_on, pars.EFTCAMB.EFTCAMB_back_turn_on)}
for alast in (0.05, 1.0):
    eta = np.array([tst * 1.0001, float(res.conformal_time(1 / alast - 1))])
    ev = res.get_time_evolution([k, k, k], eta, ["delta_photon"], lAccuracyBoost=1)
    out[str(alast)] = {"tau0": float(eta[0]), "dg": [float(v) for v in ev[:, 0, 0]]}
print("JSON" + json.dumps(out), flush=True)
"""


def child(model: str, ic: str) -> dict:
    tmp = tempfile.mkdtemp(dir=OUT)
    ini = os.path.join(tmp, "w.ini")
    with open(ini, "w") as f:
        f.write(f"DEFAULT({model if model.endswith('.ini') else str(harness.PARAMS_DIR) + '/' + model + '.ini'})\noutput_root = {tmp}/r\nfeedback_level = 0\n"
                f"EFTCAMB_write_background = F\nEFT_IC_type = {ic}\n")
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    try:
        with harness.slot(1):
            p = subprocess.run([sys.executable, "-c", CHILD, TREE, ini, str(K)], capture_output=True, text=True,
                               cwd=tmp, env=env, preexec_fn=lambda: os.nice(10), timeout=600)
    except subprocess.TimeoutExpired:
        return {"error": "timeout"}
    for line in p.stdout.splitlines():
        if line.startswith("JSON"):
            return json.loads(line[4:])
    return {"error": (p.stdout + p.stderr).strip().splitlines()[-1][:120] if (p.stdout + p.stderr).strip() else "?"}


def omega_ref(model: str, a_eval: float) -> float | None:
    f = harness.AUDIT_ROOT / "runs" / "A0" / f"{model}.npz"
    if model.endswith(".ini"):
        f = OUT / (os.path.basename(model)[:-4] + "_ic1_ref.npz")
        if not f.exists():
            harness.run_model("A0", model, f, threads=1)
    if not f.exists():
        return None
    d = harness.load(f)
    if "eft_EFTOmegaV" not in d:
        return None
    return float(np.interp(np.log(a_eval), np.log(d["eft_a_grid"]), d["eft_EFTOmegaV"]))


def one(item):
    model, ic = item
    r = child(model, ic)
    row = {"model": model, "ic": ic}
    if "error" in r:
        row["error"] = r["error"]
        return row
    x2 = (K * r["1.0"]["tau0"]) ** 2
    d05, d1 = r["0.05"]["dg"], r["1.0"]["dg"]
    row["mode1_finite"] = bool(np.isfinite(d05[0]) and np.isfinite(d1[0]))
    row["mode2_finite"] = bool(np.isfinite(d05[1]) and np.isfinite(d1[1]))
    row["stale_dep"] = abs(d05[1] / d1[1] - 1) if (d1[1] != 0 and np.isfinite(d05[1] / d1[1])) else None
    row["mode1_vs_mode2"] = abs(d1[0] / d1[1] - 1) if (d1[1] != 0 and np.isfinite(d1[0] / d1[1])) else None
    row["dg_over_x2"] = [d05[1] / x2, d1[1] / x2]
    if ic == "2":
        rr = [3 * abs(v) / (5 * x2) for v in (d05[1], d1[1])]
        row["Omega_used"] = [(5 * q - 1) / (1 - 3 * q) for q in rr]
        row["Omega_ref"] = omega_ref(model, r["a_eval"])
    return row


models = sys.argv[1:]
if models == ["all"]:
    models = [p.stem for p in harness.model_inis()]
items = [(m, ic) for m in models for ic in ("2", "3")]
with ThreadPoolExecutor(4) as ex:
    rows = list(ex.map(one, items))
for row in rows:
    print(json.dumps(row, default=lambda v: round(v, 6) if isinstance(v, float) else str(v)))
tag = "all" if len(models) > 10 else "_".join(os.path.basename(m) for m in models)[:60]
(OUT / f"sweep_{tag}.json").write_text(json.dumps(rows, indent=1, default=str))

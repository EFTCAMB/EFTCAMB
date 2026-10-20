"""Verifier check (scg-higher-derivatives-wrong): find Scaling Cubic Galileon parameter points where the positivity /
Minkowski-limit prior verdict is decided by posbound2 / tachpos (i.e. ghostpos = D1 < 0, so posbound1 > 0), and compare
the verdict of A0 (code under audit) with P_scg3 (reporter's fix) there.

Part 1: for each parameter point, posbound1/posbound2/tachpos/ghostpos recomputed from get_eft_functions output
(SCG specialisation of compute_positivity_bounds, same as scg_derivs_priors.py) on a = logspace(-6, 0, 600).
Part 2: end-to-end CAMB runs (harness.run_model) with EFT_positivity_bounds = T or EFT_minkowski_limit = T for the
points listed in E2E.
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness  # noqa: E402
from scg_derivs_priors import WORKER, bounds  # noqa: E402

SCR = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREES = {"A0": f"{SCR}/trees/A0", "FIX": f"{SCR}/work/t4b-derivs/P_scg3"}
OUTD = Path(f"{SCR}/work/verify/scg_derivs/prior_scan")
INI = Path(__file__).resolve().parents[3] / "parameters" / "5_scg.ini"
POINTS = {
    "A-0.3_default": {},
    "A+0.003": {"Scaling_Cubic_A": "0.003"},
    "A+0.01_lam10": {"Scaling_Cubic_A": "0.01", "Scaling_Cubic_lambda": "10"},
    "A-0.005": {"Scaling_Cubic_A": "-0.005"},
    "A-0.01": {"Scaling_Cubic_A": "-0.01"},
    "A+0.05_lam5": {"Scaling_Cubic_A": "0.05", "Scaling_Cubic_lambda": "5"},
}
E2E = sys.argv[1].split(",") if len(sys.argv) > 1 else []


def eft_grid(tree, over):
    tmp = tempfile.mkdtemp(prefix="ps_", dir=OUTD)
    ini = os.path.join(tmp, "scg.ini")
    extra = "".join(f"{k} = {v}\n" for k, v in over.items())
    Path(ini).write_text(f"DEFAULT({INI})\noutput_root = {tmp}/out\nEFTCAMB_write_background = F\nfeedback_level = 0\n"
                         f"EFTCAMB_feedback_level = 0\n{extra}")
    w = os.path.join(tmp, "w.py")
    Path(w).write_text(WORKER)
    with harness.slot(1):
        p = subprocess.run([sys.executable, w, tree, ini], capture_output=True, text=True, cwd=tmp,
                           env=dict(os.environ, OMP_NUM_THREADS="1"), preexec_fn=lambda: os.nice(10), timeout=900)
    line = [x for x in p.stdout.splitlines() if x.startswith("JSON")]
    if not line:
        return {"error": (p.stdout[-600:] + p.stderr[-600:])}
    return {k: np.array(v) if isinstance(v, list) else v for k, v in json.loads(line[-1][4:]).items()}


def verdicts(b, a):
    out = {}
    pos_viol = (b["posbound1"] < 0) | (b["posbound2"] < 0)
    mink_viol = (-b["ghostpos"] <= 0) | (-b["tachpos"] < 0)
    for name, m in (("positivity", pos_viol), ("minkowski", mink_viol), ("posbound1<0", b["posbound1"] < 0),
                    ("posbound2<0", b["posbound2"] < 0), ("ghostpos>=0", b["ghostpos"] >= 0),
                    ("tachpos>0", b["tachpos"] > 0)):
        out[name] = [int(m.sum()), float(a[m].min()) if m.any() else None]
    return out


def main():
    OUTD.mkdir(parents=True, exist_ok=True)
    rep = {}
    for name, over in POINTS.items():
        rep[name] = {"overrides": over}
        g = {t: eft_grid(tr, over) for t, tr in TREES.items()}
        if any("error" in g[t] for t in g):
            rep[name]["error"] = {t: g[t].get("error", "")[-300:] for t in g}
            print(f"{name:16s} background/run failed: {rep[name]['error']['A0'][-160:]!r}")
            continue
        a = g["A0"]["a"]
        b = {t: bounds(g[t]) for t in g}
        for t in g:
            rep[name][t] = verdicts(b[t], a)
        rep[name]["max_rel_diff_posbound2"] = float(np.max(np.abs(b["A0"]["posbound2"] - b["FIX"]["posbound2"])
                                                           / np.abs(b["FIX"]["posbound2"])))
        rep[name]["max_rel_diff_tachpos"] = float(np.max(np.abs(b["A0"]["tachpos"] - b["FIX"]["tachpos"])
                                                         / np.abs(b["FIX"]["tachpos"])))
        print(f"{name:16s} D1<0 at {int((b['A0']['ghostpos'] < 0).sum())}/{a.size};  "
              f"posbound2 rel.diff {rep[name]['max_rel_diff_posbound2']:.1e}, tachpos rel.diff "
              f"{rep[name]['max_rel_diff_tachpos']:.1e}")
        for t in g:
            v = rep[name][t]
            print(f"    {t:4s} violations [n, first a]: positivity {v['positivity']}  minkowski {v['minkowski']}  "
                  f"(pb1<0 {v['posbound1<0']}, pb2<0 {v['posbound2<0']}, D1>=0 {v['ghostpos>=0']}, "
                  f"tach>0 {v['tachpos>0']})")
    for name in E2E:
        over = POINTS[name]
        for t, tree in TREES.items():
            for prior in ("EFT_positivity_bounds", "EFT_minkowski_limit"):
                out = OUTD / f"e2e_{name}_{t}_{prior}.npz"
                st = harness.run_model(tree, INI, out, overrides={**over, prior: "T", "EFTCAMB_feedback_level": "2"},
                                       threads=1, timeout=900)
                log = Path(st["log"]).read_text(errors="replace")
                fl = {s: log.count(s) for s in ("Positivity bounds are violated",
                                                "Healthy Minkowski limit does not exist")}
                rep.setdefault("e2e", {})[f"{name}:{t}:{prior}"] = {"status": st["status"], "flags": fl}
                print(f"e2e {name:16s} {t:4s} {prior:22s} status={st['status']:12s} {fl}")
    (OUTD / f"scan{'_e2e' if E2E else ''}.json").write_text(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()

"""Verifier check (scg-higher-derivatives-wrong): does the wrong EFTLambdadotdot / EFTGamma2PPP change the outcome of the
optional EFT_positivity_bounds / EFT_minkowski_limit priors for the default Scaling Cubic Galileon (5_scg.ini)?

Part 1: end-to-end CAMB runs on A0 (code under audit) and P_scg3 (reporter's patched tree), priors switched on one at a
time, through harness.run_model.
Part 2: posbound1/posbound2/tachpos recomputed in Python from get_eft_functions output on a dense a grid, with the SCG
specialisation (Omega = Gamma3 = 0) of compute_positivity_bounds (06_abstract_EFTCAMB_model.f90), for both trees.
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

SCR = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREES = {"A0": f"{SCR}/trees/A0", "FIX": f"{SCR}/work/t4b-derivs/P_scg3"}
OUTD = Path(f"{SCR}/work/verify/scg_derivs/priors")
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
INI = PARAMS / "5_scg.ini"

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_background(cp)
a = np.logspace(-6, 0, 600)
_, arr = cp.EFTCAMB.get_eft_functions(res, a)
out = {"H0": cp.H0, "a": a.tolist()}
for f in ["adotoa", "Hdot", "Hdotdot", "EFTc", "EFTcdot", "EFTLambdadot", "EFTLambdadotdot", "EFTGamma1V",
          "EFTGamma2V", "EFTGamma2P", "EFTGamma2PP", "EFTGamma2PPP", "posbound1", "posbound2", "tachpos"]:
    out[f] = [float(v) for v in np.asarray(arr[f]).ravel()]
print("JSON" + json.dumps(out))
"""


def eft_grid(tree):
    tmp = tempfile.mkdtemp(prefix="pb_", dir=OUTD)
    ini = os.path.join(tmp, "scg.ini")
    Path(ini).write_text(f"DEFAULT({INI})\noutput_root = {tmp}/out\nEFTCAMB_write_background = F\nfeedback_level = 0\n"
                         "EFTCAMB_feedback_level = 0\n")
    w = os.path.join(tmp, "w.py")
    Path(w).write_text(WORKER)
    with harness.slot(1):
        p = subprocess.run([sys.executable, w, tree, ini], capture_output=True, text=True, cwd=tmp,
                           env=dict(os.environ, OMP_NUM_THREADS="1"), preexec_fn=lambda: os.nice(10), timeout=900)
    line = [x for x in p.stdout.splitlines() if x.startswith("JSON")]
    if not line:
        raise RuntimeError(p.stdout[-2000:] + p.stderr[-2000:])
    return {k: np.array(v) if isinstance(v, list) else v for k, v in json.loads(line[-1][4:]).items()}


def bounds(d):
    """SCG (Omega = Gamma3 = 0) specialisation of EFTCAMBModelPositivityBounds."""
    a, h0 = d["a"], d["H0"] / 2.99792458e8 * 1000
    ad, Hd, Hdd = d["adotoa"], d["Hdot"], d["Hdotdot"]
    g1, g2, g2p, g2pp, g2ppp = d["EFTGamma1V"], d["EFTGamma2V"], d["EFTGamma2P"], d["EFTGamma2PP"], d["EFTGamma2PPP"]
    a2 = a * a
    num = 2 * a * h0 * g2
    D1 = -2 * d["EFTc"] + 3 * a * ad * h0 * g2
    N2 = 4 * a2 * h0**2 * g1 + a * h0 * ad * (-3 * g2 + a * g2p)
    N3 = -2 * d["EFTcdot"] + 3 * a * h0 * (g2 * Hd + ad**2 * (-g2 + a * g2p))
    N4 = -2 * d["EFTLambdadot"] * ad + 2 * d["EFTLambdadotdot"] - a2 * ad * h0 * Hd * (g2p - 3 * a * g2pp)
    N5 = a2 * h0 * Hdd * g2p + a2**2 * h0 * ad**3 * g2ppp
    pb1 = -(4 * a2 * h0**2 * g2**2) / D1 / 8
    pb2 = 0.5 * (N2 + num * -N3 / D1 + num**2 * (-N4 - N5) / D1**2 / 8)
    Q1 = (-ad * d["EFTLambdadot"] + d["EFTLambdadotdot"]) * a2
    Q2 = h0 * (Hdd - ad * Hd) * a2 * g2p / 2
    Q3 = 3 * h0 * ad * Hd * a**3 * g2pp / 2
    Q4 = h0 * ad**3 * a2**2 * g2ppp / 2
    tach = -(Q1 + Q2 + Q3 + Q4) / D1
    return {"ghostpos": D1, "posbound1": pb1, "posbound2": pb2, "tachpos": tach}


def main():
    OUTD.mkdir(parents=True, exist_ok=True)
    rep = {"runs": {}, "grid": {}}
    # part 1: end-to-end runs with the priors on
    for tag, tree in TREES.items():
        for prior in ("EFT_positivity_bounds", "EFT_minkowski_limit"):
            out = OUTD / f"{tag}_{prior}.npz"
            st = harness.run_model(tree, INI, out, overrides={prior: "T", "EFTCAMB_feedback_level": "2"}, threads=1,
                                   timeout=900)
            log = Path(st["log"]).read_text(errors="replace") if st.get("log") else ""
            flags = {s: log.count(s) for s in ("Positivity bounds are violated", "Healthy Minkowski limit does not exist")}
            meta = {}
            if out.exists():
                meta = harness.load(out)["meta"]
            rep["runs"][f"{tag}:{prior}"] = {"status": st["status"], "log_flags": flags,
                                             "errors": meta.get("errors"), "tail": "\n".join(st.get("log_tail", []))[-400:]}
            print(f"{tag:4s} {prior:22s} status={st['status']:12s} flags={flags} errors={str(meta.get('errors'))[:200]}")
    # part 2: bound values on a dense grid
    g = {tag: eft_grid(tree) for tag, tree in TREES.items()}
    b = {tag: bounds(g[tag]) for tag in g}
    a = g["A0"]["a"]
    for key in ("ghostpos", "posbound1", "posbound2", "tachpos"):
        sa, sf = np.sign(b["A0"][key]), np.sign(b["FIX"][key])
        dif = np.where(sa != sf)[0]
        relmax = float(np.max(np.abs(b["A0"][key] - b["FIX"][key]) / np.maximum(np.abs(b["FIX"][key]), 1e-300)))
        crit = {"posbound1": "<0", "posbound2": "<0", "tachpos": "-x<0", "ghostpos": "-x<=0"}[key]
        viol = {}
        for tag in ("A0", "FIX"):
            x = b[tag][key]
            m = {"posbound1": x < 0, "posbound2": x < 0, "tachpos": -x < 0, "ghostpos": -x <= 0}[key]
            viol[tag] = [float(a[m].min()), float(a[m].max()), int(m.sum())] if m.any() else None
        rep["grid"][key] = {"sign_differs_n": int(dif.size), "a_sign_differs": [float(a[dif].min()), float(a[dif].max())]
                            if dif.size else None, "max_rel_diff": relmax, "violated_A0": viol["A0"],
                            "violated_FIX": viol["FIX"], "criterion": crit}
        print(f"grid {key:10s} max|A0-FIX|/|FIX|={relmax:.2e}  sign differs at {dif.size}/{a.size} points"
              f" {rep['grid'][key]['a_sign_differs']}; violated(A0)={viol['A0']} violated(FIX)={viol['FIX']}")
        for ai in (1e-3, 1e-2, 0.1, 0.5, 1.0):
            i = int(np.argmin(np.abs(np.log(a / ai))))
            print(f"     a={a[i]:.3g}: A0={b['A0'][key][i]:+.4e} FIX={b['FIX'][key][i]:+.4e}")
    # sanity: python-recomputed tachpos vs Fortran field if filled
    rep["fortran_fields_filled"] = {k: bool(np.any(g["A0"][k] != 0)) for k in ("posbound1", "posbound2", "tachpos")}
    print("fortran posbound fields filled by get_eft_functions:", rep["fortran_fields_filled"])
    (OUTD / "priors.json").write_text(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()

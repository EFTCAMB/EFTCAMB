"""Accuracy of the h=0.08 EFTLambda-value stencil used in the P_val diagnostic build vs node-based / EFTLambdadot-based
second derivatives, for 5_quint_exponential (where do the Minkowski fail counts differ?).  Usage: python minktach_stencil.py"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

os.environ["OMP_NUM_THREADS"] = "1"
ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
sys.path.insert(0, f"{ROOT}/trees/A0")
PAR = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"
with harness.slot(1):
    import numpy as np
    import camb
    cp = camb.read_ini(f"{PAR}/5_quint_exponential.ini", no_validate=True)
    cp.EFTCAMB.EFTCAMB_skip_stability = True
    cp.EFTCAMB.EFTCAMB_feedback_level = 0
    cp.EFTCAMB.EFTCAMB_write_background = False
    res = camb.get_background(cp)
    eft = cp.EFTCAMB if False else res.Params.EFTCAMB
    a = np.logspace(-3, -0.05, 2000)
    h = 0.08
    get = lambda aa: eft.get_eft_functions(res, aa)[1]
    c0, cp_, cm = get(a), get(a * np.exp(h)), get(a * np.exp(-h))
    L = lambda c, aa: np.array(c["EFTLambda"]) / aa**2
    l0, lp, lm = L(c0, a), L(cp_, a * np.exp(h)), L(cm, a * np.exp(-h))
    H, Hd, Ld, c = (np.array(c0[f]) for f in ("adotoa", "Hdot", "EFTLambdadot", "EFTc"))
    Ldd_h = a**2 * (Hd * (lp - lm) / (2 * h) + H**2 * (lp - 2 * l0 + lm) / h**2)
    # EFTLambdadot-based with the same h
    Ldp, Ldm = np.array(cp_["EFTLambdadot"]), np.array(cm["EFTLambdadot"])
    Ldd_yp = H * (Ldp - Ldm) / (2 * h) - 2 * H * Ld
    t = lambda Ldd: (Ldd - H * Ld) * a**2 / (2 * c)
    for name, x in (("value h=0.08", Ldd_h), ("EFTLambdadot h=0.08", Ldd_yp)):
        print(f"{name:22s} fails {(t(x) > 0).sum()}/{len(a)}")
    rel = np.abs(Ldd_h - Ldd_yp) / np.abs(Ldd_yp)
    print("median |value-yp|/|yp| =", float(np.median(rel)), " 90th pct =", float(np.percentile(rel, 90)))
    print("ratio (Ldd - H Ld)/Ldd at a = 0.001, 0.01, 0.1, 0.5:",
          [float(((Ldd_yp - H * Ld) / Ldd_yp)[np.argmin(abs(a - q))]) for q in (1e-3, 1e-2, 0.1, 0.5)])
    dis = (t(Ldd_h) > 0) != (t(Ldd_yp) > 0)
    print("sign disagreements:", int(dis.sum()), "a range", (float(a[dis].min()), float(a[dis].max())) if dis.any() else None)

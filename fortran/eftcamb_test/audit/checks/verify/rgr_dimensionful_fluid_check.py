"""Verifier: compare EFTCAMB 4_Des_mc5e_CPL_1 at GR thresholds 1e-8 / 1e-7 (A0) against plain CAMB v2 (B0) with
a CPL fluid (cs2=1, perturbations on), and scan the threshold to locate the flip point vs the dimensionful
EFT functions. Independent of the repro: different code (upstream DE fluid) and a threshold scan."""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "rgr_dimensionful"
INI = harness.PARAMS_DIR / "4_Des_mc5e_CPL_1.ini"

WORKER = r"""
import sys, json, numpy as np
cfg = json.loads(sys.argv[1])
sys.path.insert(0, cfg["tree"])
import camb
pars = camb.read_ini(cfg["ini"])
pars.set_dark_energy(w=-0.9, wa=-0.1, cs2=1.0, dark_energy_model=cfg["de"])
r = camb.get_results(pars)
cl = r.get_total_cls(2500, CMB_unit="muK", raw_cl=True)
pp = r.get_lens_potential_cls(2500, raw_cl=True)
np.savez(cfg["out"], cl_total=cl, cl_lens_potential=pp, sigma8=r.get_sigma8())
print("camb from", camb.__file__)
"""


def run_b0(de):
    out = OUT / f"B0_{de}.npz"
    wf = OUT / "b0_worker.py"
    wf.write_text(WORKER)
    cfg = {"tree": str(harness.tree_path("B0")), "ini": str(INI), "out": str(out), "de": de}
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, str(wf), json.dumps(cfg)], cwd=OUT, env=env,
                           capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(p.stderr[-2000:])
    print(p.stdout.strip())
    return harness.load(out)


def run_a0(thr):
    out = OUT / f"A0_thr{thr}.npz"
    st = harness.run_model("A0", INI, out, overrides={"EFTCAMB_GR_threshold": thr}, threads=1, timeout=600)
    assert st["status"] == "ok", st
    rgr = re.findall(r"Return to GR time:\s*(\S+)", out.with_suffix(".log").read_text(errors="replace"))
    return harness.load(out), float(rgr[0])


def stats(d, ref):
    tt = d["cl_total"][2:2501, 0] / ref["cl_total"][2:2501, 0] - 1
    pp = d["cl_lens_potential"][2:2501, 0] / ref["cl_lens_potential"][2:2501, 0] - 1
    return (f"dsig8={float(np.ravel(d['sigma8'])[0]) / float(np.ravel(ref['sigma8'])[0]) - 1:+.2e} "
            f"dTT2={tt[0]:+.2e} maxdTT30={np.abs(tt[28:]).max():.2e} maxdPP={np.abs(pp).max():.2e}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fluid = run_b0("fluid")
    ppf = run_b0("ppf")
    res = {}
    for thr in ("1e-8", "5.5e-9", "6e-9", "2e-8", "5e-8", "1e-7"):
        res[thr] = run_a0(thr)
        print(f"A0 thr={thr}: RGR={res[thr][1]:.4g} sigma8={float(np.ravel(res[thr][0]['sigma8'])[0]):.5f}")
    d0 = res["1e-8"][0]
    a = d0["eft_a_grid"]
    for f in ("EFTc", "EFTLambda", "EFTcdot", "EFTLambdadot", "EFTOmegaV", "EFTOmegaP"):
        k = "eft_" + f
        if k in d0:
            print(f"  max |{f}| raw over a-grid: {np.nanmax(np.abs(d0[k])):.3e}; /a^2: {np.nanmax(np.abs(d0[k] / a**2)):.3e}")
    print("B0 fluid sigma8 =", float(np.ravel(fluid["sigma8"])[0]), " B0 ppf sigma8 =", float(np.ravel(ppf["sigma8"])[0]))
    print("ppf vs fluid:", stats(ppf, fluid))
    for thr, (d, r) in res.items():
        print(f"A0 thr={thr} (RGR={r:.3g}) vs B0 fluid:", stats(d, fluid))


if __name__ == "__main__":
    main()

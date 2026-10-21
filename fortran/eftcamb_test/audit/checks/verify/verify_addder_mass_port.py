"""Independent check of finding additional-derivs-wrong-cache.

Method (independent of the t2d-static repro and of the FIX1 build):
  1. For every shipped ini with EFT_mass_stability = T, run the model on a tree (default A0) with the mass check off,
     then on a grid of a in (EFTCAMB_pert_turn_on, 1] call the public EFTTestStability with EFT_mass_stability = T
     (this runs compute_additional_derivs + compute_mass_factors) and keep the whole timestep cache.
  2. Get the model's own analytic derivatives with get_eft_functions (a different Fortran code path that never calls
     compute_additional_derivs) on the same grid and at a*(1 +- eps); check analytic second derivatives against a centred
     finite difference of the analytic first derivatives (validates the reference and the dot convention).
  3. Port compute_mass_factors to Python by mechanically translating the Fortran source of
     06_abstract_EFTCAMB_model.f90 (no hand transcription), validate the port on the code's own cache (must reproduce the
     code's EFT_mu1/EFT_mu2), then re-evaluate mu1/mu2 with EFTcdotdot/EFTGamma{1,2,3}PP replaced by the analytic values.
  4. Report the verdict of the mass criterion (mu < -(rate*H/a)^2 for a > turn-on) with code vs corrected inputs.

Usage: python verify_addder_mass_port.py OUTDIR [TREE] [TREE2 ...]
Each tree is processed; for trees after the first only the code values are reported (e.g. FIX1 for comparison).
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REPO = HERE.parents[4]
SCR = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
PAR = REPO / "fortran/eftcamb_test/parameters"
SRC = REPO / "fortran/eftcamb/06_abstract_EFTCAMB_model.f90"
REPLACED = ["EFTcdotdot", "EFTGamma1PP", "EFTGamma2PP", "EFTGamma3PP"]

CHILD = r"""
import sys, json, ctypes, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
from camb import eftcamb
from camb.baseconfig import camblib, CAMBError
ini, out = sys.argv[2], sys.argv[3]
res_out = {}
pars = camb.read_ini(ini, no_validate=True)
try:
    camb.get_background(pars); res_out["verdict"] = "stable"
except CAMBError as e:
    res_out["verdict"] = str(e).strip().splitlines()[-1]
pars = camb.read_ini(ini, no_validate=True)
pars.EFTCAMB.EFT_mass_stability = False
pars.EFTCAMB.EFTCAMB_feedback_level = 0
res = camb.get_background(pars)
P = res.Params
ton = P.EFTCAMB.EFTCAMB_pert_turn_on
rate = P.EFTCAMB.EFT_mass_stability_rate
a = np.unique(np.concatenate([np.linspace(ton, 1.0, 400)[1:], ton + (1 - ton) * np.logspace(-6, -1, 100)]))
eps = 1e-4
P.EFTCAMB.EFT_mass_stability = True
f = camblib.__eftcamb_stability_MOD_eftteststability
f.restype = ctypes.c_int
arr = (eftcamb.EFTCAMB_timestep_cache * a.size)()
ok = np.zeros(a.size, dtype=bool)
for i, ai in enumerate(a):
    kmax = ctypes.c_double(10.0)
    ok[i] = bool(f(ctypes.byref(ctypes.c_double(ai)), ctypes.byref(kmax), ctypes.byref(P.EFTCAMB),
                   ctypes.byref(P.EFTCAMB_parameter_cache), ctypes.byref(arr[i])))
used = np.ctypeslib.as_array(arr)
P.EFTCAMB.EFT_mass_stability = False
_, ana = P.EFTCAMB.get_eft_functions(res, a)
_, anap = P.EFTCAMB.get_eft_functions(res, a * (1 + eps))
_, anam = P.EFTCAMB.get_eft_functions(res, a * (1 - eps))
np.savez(out, a=a, ok=ok, used=used, ana=ana, anap=anap, anam=anam, eps=eps, ton=ton, rate=rate,
         h0_Mpc=P.EFTCAMB_parameter_cache.h0_Mpc)
print("JSON" + json.dumps(res_out))
"""


def fortran_mass_factors():
    """Translate the body of EFTCAMBModelComputeMassFactors into Python source (vectorised over a)."""
    lines = SRC.read_text().splitlines()
    i0 = next(i for i, ln in enumerate(lines) if "subroutine EFTCAMBModelComputeMassFactors" in ln)
    i1 = next(i for i, ln in enumerate(lines) if "end subroutine EFTCAMBModelComputeMassFactors" in ln)
    body = lines[i0:i1]
    start = next(i for i, ln in enumerate(body) if "real(dl) :: EFTcP, EFTcPP" in ln) + 1
    stmts, cur = [], ""
    for ln in body[start:]:
        s = ln.split("!")[0].rstrip()
        if not s.strip():
            continue
        if s.endswith("&"):
            cur += s[:-1].strip() + " "
            continue
        cur += s.strip()
        stmts.append(cur)
        cur = ""
    py = []
    for st in stmts:
        st = re.sub(r"(\d+\.?\d*)_dl", r"\1", st)
        st = st.replace("eft_par_cache%h0_Mpc", "h0").replace("eft_cache%EFT_mu", "OUT_mu")
        st = re.sub(r"eft_cache%(\w+)", r"c['\1']", st)
        st = re.sub(r"\bsqrt\(", "np.sqrt(", st)
        assert "%" not in st, st
        py.append(st)
    return "\n".join(py), len(stmts)


def mass_factors(code_src, c, a, h0):
    ns = {"np": np, "c": c, "a": a, "h0": h0}
    with np.errstate(all="ignore"):
        exec(code_src, ns)  # noqa: S102
    return ns["OUT_mu1"], ns["OUT_mu2"]


def run_tree(tree, ini, outdir):
    out = outdir / f"{Path(tree).name}_{Path(ini).stem}.npz"
    w = outdir / f"w_{Path(ini).stem}.ini"
    w.write_text(f"DEFAULT({ini})\nfeedback_level = 0\nEFTCAMB_feedback_level = 0\nEFTCAMB_write_background = F\n"
                 f"output_root = {outdir / Path(ini).stem}\n")
    with harness.slot(1):
        r = subprocess.run([sys.executable, "-c", CHILD, str(tree), str(w), str(out)], capture_output=True, text=True,
                           cwd=outdir, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900,
                           preexec_fn=lambda: os.nice(10))
    js = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
    if not js:
        raise RuntimeError(f"{tree} {ini}: rc={r.returncode}\n{r.stdout[-800:]}\n{r.stderr[-800:]}")
    return json.loads(js[-1][4:]), dict(np.load(out, allow_pickle=False))


def todict(rec):
    return {k: np.asarray(rec[k], dtype=float) for k in rec.dtype.names}


def main():
    outdir = Path(sys.argv[1])
    outdir.mkdir(parents=True, exist_ok=True)
    trees = [Path(t) for t in sys.argv[2:]] or [SCR / "trees/A0"]
    inis = sorted(p for p in PAR.glob("*.ini")
                  if re.search(r"^\s*EFT_mass_stability\s*=\s*T", p.read_text(), re.M | re.I))
    src, nst = fortran_mass_factors()
    print(f"translated {nst} Fortran statements of EFTCAMBModelComputeMassFactors; {len(inis)} inis with mass check")
    summary = []
    for ini in inis:
        row = {"model": ini.stem}
        for it, tree in enumerate(trees):
            verdict, d = run_tree(tree, ini, outdir)
            a, h0, ton, rate = d["a"], float(d["h0_Mpc"]), float(d["ton"]), float(d["rate"])
            used, ana = todict(d["used"]), todict(d["ana"])
            anap, anam = todict(d["anap"]), todict(d["anam"])
            H = used["adotoa"]
            thr = -(rate * H / a) ** 2
            mu1c, mu2c = used["EFT_mu1"], used["EFT_mu2"]
            tag = Path(tree).name
            row[f"verdict_{tag}"] = verdict["verdict"]
            row[f"code_unstable_grid_{tag}"] = bool(np.any(((mu1c < thr) | (mu2c < thr)) & (a > ton)))
            row[f"min_mu1_over_H2_{tag}"] = float(np.nanmin(mu1c * (a / H) ** 2))
            if it > 0:
                continue
            # port validation on the code's own cache
            p1, p2 = mass_factors(src, used, a, h0)
            row["port_relerr_mu1"] = float(np.nanmax(np.abs(p1 - mu1c) / np.maximum(np.abs(mu1c), 1e-300)))
            row["port_relerr_mu2"] = float(np.nanmax(np.abs(p2 - mu2c) / np.maximum(np.abs(mu2c), 1e-300)))
            # analytic reference check: second derivative vs centred FD of analytic first derivative
            hstep = 2 * a * float(d["eps"])
            fd = {
                "EFTGamma1PP": (anap["EFTGamma1P"] - anam["EFTGamma1P"]) / hstep,
                "EFTGamma2PP": (anap["EFTGamma2P"] - anam["EFTGamma2P"]) / hstep,
                "EFTGamma3PP": (anap["EFTGamma3P"] - anam["EFTGamma3P"]) / hstep,
                "EFTcdotdot": a * ana["adotoa"] * (anap["EFTcdot"] - anam["EFTcdot"]) / hstep
                - 2 * ana["adotoa"] * ana["EFTcdot"],
            }
            for k in REPLACED:
                scale = max(np.max(np.abs(ana[k])), 1e-300)
                row[f"ana_vs_fd_{k}"] = float(np.max(np.abs(ana[k] - fd[k])) / scale)
                row[f"used_vs_ana_{k}"] = float(np.max(np.abs(used[k] - ana[k])) / scale) if scale > 1e-300 else 0.0
            row["used_cdotdot_eq_m2Hcdot"] = float(np.max(np.abs(used["EFTcdotdot"] + 2 * H * used["EFTcdot"]))
                                                   / max(np.max(np.abs(used["EFTcdotdot"])), 1e-300))
            corr = dict(used)
            for k in REPLACED:
                corr[k] = ana[k]
            q1, q2 = mass_factors(src, corr, a, h0)
            row["port_corrected_unstable_grid"] = bool(np.any(((q1 < thr) | (q2 < thr)) & (a > ton)))
            row["port_corrected_min_mu1_over_H2"] = float(np.nanmin(q1 * (a / H) ** 2))
            row["port_corrected_min_mu2_over_H2"] = float(np.nanmin(q2 * (a / H) ** 2))
            np.savez(outdir / f"port_{ini.stem}.npz", a=a, mu1_code=mu1c, mu2_code=mu2c, mu1_corr=q1, mu2_corr=q2,
                     H=H)
        summary.append(row)
        print(json.dumps(row))
    (outdir / "summary.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

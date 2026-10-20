"""Independent check of finding ini-theta-h0-solved-without-eft.

Claim: with thetastar / cosmomc_theta in an ini, CAMB_ReadParams solves H0 with P%EFTCAMB still defaulted (EFTflag=0),
so for EFT models with non-LCDM background the H0 does not reproduce the requested theta.

Independent method (no set_cosmology, no derived thetastar of the theta-ini run itself):
  A. read_ini(ini with thetastar)                       -> H0_ini
  B. read_ini(ini with hubble=H0_ini, EFT on)            -> theta_EFT(H0_ini)   (expected != requested if bug)
  C. read_ini(ini with hubble=H0_ini, EFTflag=0)         -> theta_GR(H0_ini)    (== requested if solve used EFTflag=0)
  D. own scipy brentq on hubble with EFT on (hubble-ini path only) -> H0_true
  E. CAMB_RunFromIni path (camb.run_ini, subprocess, FeedbackLevel=1): parse printed H0 and 100*theta
  F. cosmomc_theta variant (read_ini) for selected models: results.cosmomc_theta() vs requested
Run: /Users/marco/Software/miniforge3/bin/python ini_theta_indep.py [models...]
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

TREE = harness.tree_path("A0") if hasattr(harness, "tree_path") else None
TREE = str(TREE) if TREE else (
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0")
OUT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit/work/verify/ini_theta")
OUT.mkdir(parents=True, exist_ok=True)
PARAMS = HERE.parents[2] / "parameters"
TREE = os.environ.get("AUDIT_TREE", TREE)
sys.path.insert(0, TREE)
import camb  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

THETA = 0.0104110
CTHETA = 0.0104092


def flatten(path: Path) -> dict:
    keys, defaults = {}, []
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.upper().startswith("DEFAULT("):
            defaults.append(flatten(path.parent / line[8:line.rindex(")")].strip()))
        elif "=" in line:
            k, v = line.split("=", 1)
            keys[k.strip()] = v.strip()
    merged = {}
    for d in defaults:
        merged.update(d)
    return merged | keys


def write_ini(name: str, keys: dict) -> Path:
    p = OUT / f"{name}.ini"
    p.write_text("\n".join(f"{k} = {v}" for k, v in keys.items()) + "\n")
    return p


CHEAP = {"get_scalar_cls": "F", "get_tensor_cls": "F", "get_transfer": "F", "num_redshiftwindows": "0",
         "feedback_level": "0", "EFTCAMB_write_background": "F"}


def theta_from_hubble(model: str, base: dict, H0: float, eftflag0: bool = False, kind: str = "thetastar") -> float:
    keys = dict(base) | CHEAP | {"hubble": repr(float(H0))}
    if eftflag0:
        keys["EFTflag"] = "0"
    cp = camb.read_ini(str(write_ini(f"{model}_hub", keys)))
    if eftflag0:
        assert cp.EFTCAMB.EFTflag == 0
    res = camb.get_background(cp)
    if kind == "thetastar":
        return res.get_derived_params()["thetastar"] / 100
    return res.cosmomc_theta()


def run_ini_path(model: str, base: dict) -> dict:
    keys = dict(base) | {"thetastar": str(THETA), "get_scalar_cls": "T", "get_tensor_cls": "F", "get_transfer": "F",
                         "do_lensing": "F", "l_max_scalar": "300", "num_redshiftwindows": "0",
                         "EFTCAMB_write_background": "F", "output_root": str(OUT / f"runini_{model}")}
    ini = write_ini(f"{model}_runini", keys)
    code = (f"import sys; sys.path.insert(0, {TREE!r}); import camb; camb.set_feedback_level(1); "
            f"camb.run_ini({str(ini)!r})")
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=500,
                       env=dict(os.environ, OMP_NUM_THREADS="1"))
    (OUT / f"runini_{model}.log").write_text(p.stdout + "\n---stderr---\n" + p.stderr)
    h0 = re.findall(r"^\s*H0\s*=\s*([0-9.]+)", p.stdout, re.M)
    th = re.findall(r"100\*theta\s*=\s*([0-9.]+)", p.stdout)
    return {"rc": p.returncode, "H0": float(h0[-1]) if h0 else None, "theta": float(th[-1]) / 100 if th else None}


def main(models):
    results = {}
    for model in models:
        base = flatten(PARAMS / f"{model}.ini")
        base.pop("hubble")
        r = {"EFTflag": base.get("EFTflag")}
        # A: ini thetastar path (CAMB_ReadParamFile)
        cp = camb.read_ini(str(write_ini(f"{model}_th", dict(base) | CHEAP | {"thetastar": str(THETA)})))
        r["H0_ini"] = cp.H0
        # B, C
        r["theta_EFT_at_H0ini"] = theta_from_hubble(model, base, cp.H0)
        r["theta_GR_at_H0ini"] = theta_from_hubble(model, base, cp.H0, eftflag0=True)
        # D: own root find in hubble with EFT on
        lo = max(40.0, 1.1 * 100 * (float(base["ombh2"]) + float(base["omch2"]) + float(base.get("omnuh2", 0))) ** 0.5)
        try:
            r["H0_true"] = brentq(lambda h: theta_from_hubble(model, base, h) - THETA, lo, 100.0, rtol=1e-6)
        except ValueError:
            r["H0_true"] = float("nan")
            r["scan"] = {h: theta_from_hubble(model, base, h) / THETA - 1 for h in (lo, 50.0, 60.0, 70.0, 80.0, 100.0)}
            print("  no bracket; theta_EFT(H0)/req-1 scan:", r["scan"], flush=True)
        r["rel_theta_EFT_at_H0ini"] = r["theta_EFT_at_H0ini"] / THETA - 1
        r["rel_theta_GR_at_H0ini"] = r["theta_GR_at_H0ini"] / THETA - 1
        # E: CAMB_RunFromIni path
        r["run_ini"] = run_ini_path(model, base)
        if r["run_ini"]["theta"]:
            r["run_ini"]["rel_theta"] = r["run_ini"]["theta"] / THETA - 1
        print(f"{model:22s} EFTflag={r['EFTflag']} H0_ini={r['H0_ini']:.4f} H0_true={r['H0_true']:.4f} "
              f"theta_EFT(H0_ini) rel={r['rel_theta_EFT_at_H0ini']:+.3e} theta_GR(H0_ini) rel="
              f"{r['rel_theta_GR_at_H0ini']:+.3e} | run_ini: H0={r['run_ini']['H0']} "
              f"rel={r['run_ini'].get('rel_theta')}", flush=True)
        results[model] = r
    return results


def cosmomc(models):
    results = {}
    for model in models:
        base = flatten(PARAMS / f"{model}.ini")
        base.pop("hubble")
        cp = camb.read_ini(str(write_ini(f"{model}_cth", dict(base) | CHEAP | {"cosmomc_theta": str(CTHETA)})))
        got = camb.get_background(cp).cosmomc_theta()
        gr = theta_from_hubble(model, base, cp.H0, eftflag0=True, kind="cosmomc")
        results[model] = {"H0_ini": cp.H0, "rel": got / CTHETA - 1, "rel_GR_at_H0ini": gr / CTHETA - 1}
        print(f"cosmomc_theta {model:22s} H0_ini={cp.H0:.4f} rel={got / CTHETA - 1:+.3e} "
              f"GR-at-H0ini rel={gr / CTHETA - 1:+.3e}", flush=True)
    return results


if __name__ == "__main__":
    args = sys.argv[1:]
    mode = "theta"
    if args and args[0] == "--cosmomc":
        mode, args = "cosmomc", args[1:]
    with harness.slot(1):
        res = main(args) if mode == "theta" else cosmomc(args)
    tag = Path(TREE).name + "_" + "_".join(args)[:80]
    (OUT / f"indep_{mode}_{tag}.json").write_text(json.dumps(res, indent=1, default=str))

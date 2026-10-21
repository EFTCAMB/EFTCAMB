"""Adversarial verification of t3-get-evolution-transfer-only-segfault (independent of the t3-gr repro).

Differences from the reporter's repro:
  * only cp.WantCls = False is set (WantScalars is left as the ini sets it, T): CAMB_GetResults itself
    forces P%WantScalars = .false. for the transfer-only pass, i.e. the standard matter-power-only workflow;
  * EFTCAMB.get_evolution is called directly with conformal times (not via get_scale_evolution);
  * the upstream entry point results.get_time_evolution is run on the same results object;
  * other trees: B0 (plain CAMB v2, upstream get_time_evolution), B1 (pre-rebase EFTCAMB),
    an -fcheck=all build with only F-002 patched (work/t1-crash/chk/AD) to locate the faulting access;
  * a second EFT model (4_Des_F_R_LCDM_1) to rule out model specificity.
Each case runs in its own subprocess holding one harness CPU slot. Output: work/verify/getevol_indep.json
"""

import json
import os
import signal
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "getevol"
CHK_AD = harness.AUDIT_ROOT / "work" / "t1-crash" / "chk" / "AD"

WORKER = r"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
tree, ini, mode, wantcls, wantscal, wanttens, checks_dir = sys.argv[1:8]
sys.path.insert(0, tree)
import camb, numpy as np
sys.path.append(checks_dir)
assert os.path.realpath(camb.__file__).startswith(os.path.realpath(tree)), camb.__file__
wantcls_arg = wantcls
wantcls = wantcls == "T"
if ini == "GR":
    cp = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965, WantCls=wantcls,
                         WantTransfer=True, kmax=2.0, lmax=300)
else:
    cp = camb.read_ini(ini, no_validate=True)
    if hasattr(cp, "EFTCAMB") and cp.EFTCAMB.EFTflag == 0:
        import harness
        keys = harness.parse_ini(ini)
        if int(keys.get("EFTflag", "0")) != 0:
            cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)
            print("B1 fallback EFT init", flush=True)
if wantcls_arg != "keep":
    cp.WantCls = wantcls
if wantscal != "keep":
    cp.WantScalars = wantscal == "T"
if wanttens != "keep":
    cp.WantTensors = wanttens == "T"
if ini != "GR":
    if not cp.WantTensors:
        cp.InitPower.r = 0.0
    elif cp.InitPower.r == 0:
        cp.InitPower.r = 0.1
    cp.DoLensing = False
    cp.max_l = 300
    cp.WantTransfer = True
    cp.SourceWindows = []
print("before get_results: WantCls", cp.WantCls, "WantScalars", cp.WantScalars, "WantTensors", cp.WantTensors,
      "EFTflag", getattr(getattr(cp, "EFTCAMB", None), "EFTflag", None), flush=True)
r = camb.get_results(cp)
print("get_results ok; results.Params.WantScalars =", bool(r.Params.WantScalars), "WantCls =", bool(r.Params.WantCls),
      "WantTensors =", bool(r.Params.WantTensors), flush=True)
eta = np.asarray(r.conformal_time(np.array([9.0, 1.0, 0.0])))
if mode == "eft":
    names, arr = cp.EFTCAMB.get_evolution(r, [0.05], eta)
    print("RESULT eft pi =", json.dumps(np.asarray(arr["pi"])[0].tolist()), flush=True)
else:
    ev = r.get_time_evolution([0.05], eta, ["delta_cdm", "Weyl"])
    print("RESULT upstream delta_cdm,Weyl =", json.dumps(np.asarray(ev)[0].tolist()), flush=True)
"""

CASES = [
    # name, tree, ini, mode, WantCls[, WantScalars, WantTensors]  ("keep" = leave as ini/set_params)
    ("A0_RPH_MK_cls_eft", "A0", "3_RPH_MK_const_1", "eft", "T"),
    ("A0_RPH_MK_nocls_eft", "A0", "3_RPH_MK_const_1", "eft", "F"),
    ("A0_RPH_MK_nocls_eft_rerun", "A0", "3_RPH_MK_const_1", "eft", "F"),
    ("A0_DesFR_nocls_eft", "A0", "4_Des_F_R_LCDM_1", "eft", "F"),
    ("A0_RPH_MK_cls_upstream", "A0", "3_RPH_MK_const_1", "up", "T"),
    ("A0_RPH_MK_nocls_upstream", "A0", "3_RPH_MK_const_1", "up", "F"),
    ("A0_GR_nocls_upstream", "A0", "GR", "up", "F"),
    ("B0_GR_cls_upstream", "B0", "GR", "up", "T"),
    ("B0_GR_nocls_upstream", "B0", "GR", "up", "F"),
    ("B1_RPH_MK_cls_eft", "B1", "3_RPH_MK_const_1", "eft", "T"),
    ("B1_RPH_MK_nocls_eft", "B1", "3_RPH_MK_const_1", "eft", "F"),
    ("B1_RPH_MK_nocls_upstream", "B1", "3_RPH_MK_const_1", "up", "F"),
    ("chkAD_RPH_MK_cls_eft", str(CHK_AD), "3_RPH_MK_const_1", "eft", "T"),
    ("chkAD_RPH_MK_nocls_eft", str(CHK_AD), "3_RPH_MK_const_1", "eft", "F"),
    ("chkAD_RPH_MK_nocls_upstream", str(CHK_AD), "3_RPH_MK_const_1", "up", "F"),
    # explicit WantScalars = F (as in the reporter's repro)
    ("A0_RPH_MK_noscal_eft", "A0", "3_RPH_MK_const_1", "eft", "F", "F", "F"),
    ("A0_RPH_MK_noscal_upstream", "A0", "3_RPH_MK_const_1", "up", "F", "F", "F"),
    ("A0_DesFR_noscal_eft", "A0", "4_Des_F_R_LCDM_1", "eft", "F", "F", "F"),
    ("A0_GR_noscal_upstream", "A0", "GR", "up", "F", "F", "F"),
    ("B0_GR_noscal_upstream", "B0", "GR", "up", "F", "F", "F"),
    ("B1_RPH_MK_noscal_eft", "B1", "3_RPH_MK_const_1", "eft", "F", "F", "F"),
    ("B1_RPH_MK_noscal_upstream", "B1", "3_RPH_MK_const_1", "up", "F", "F", "F"),
    # ini only: get_scalar_cls = F in a wrapper ini, no Python attribute changes except output trimming
    ("A0_RPH_MK_ininoscal_eft", "A0", "INI:3_RPH_MK_const_1", "eft", "keep", "keep", "keep"),
    # tensor-only C_l run (WantCls=T, WantScalars=F, WantTensors=T)
    ("A0_RPH_MK_tensonly_eft", "A0", "3_RPH_MK_const_1", "eft", "T", "F", "T"),
    ("chkAD_RPH_MK_noscal_eft", str(CHK_AD), "3_RPH_MK_const_1", "eft", "F", "F", "F"),
    ("chkAD_RPH_MK_noscal_upstream", str(CHK_AD), "3_RPH_MK_const_1", "up", "F", "F", "F"),
    # 3_RPH_MK_const_1 trips an unrelated -fcheck bound in 02_equispaced_interpolation_linear_1D: use other models
    ("chkAD_DesFR_cls_eft", str(CHK_AD), "4_Des_F_R_LCDM_1", "eft", "T", "T", "F"),
    ("chkAD_DesFR_noscal_eft", str(CHK_AD), "4_Des_F_R_LCDM_1", "eft", "F", "F", "F"),
    ("chkAD_GR_cls_upstream", str(CHK_AD), "GR", "up", "T", "T", "F"),
    ("chkAD_GR_noscal_upstream", str(CHK_AD), "GR", "up", "F", "F", "F"),
    # set_params trips F-001 on the -fcheck tree: GR (EFTflag=0) through read_ini instead
    ("chkAD_base_noscal_upstream", str(CHK_AD), "base_params", "up", "F", "F", "F"),
    ("chkAD_base_noscal_eft", str(CHK_AD), "base_params", "eft", "F", "F", "F"),
]


def run_case(name, tree, ini, mode, wantcls, wantscal="keep", wanttens="F"):
    OUT.mkdir(parents=True, exist_ok=True)
    if ini == "GR":
        ini_path = ini
    elif ini.startswith("INI:"):
        ini_path = str(harness._wrapper_ini(harness.PARAMS_DIR / f"{ini[4:]}.ini", OUT, name,
                                            {"get_scalar_cls": "F", "get_tensor_cls": "F", "get_vector_cls": "F",
                                             "get_transfer": "T", "do_lensing": "F"}))
    else:
        ini_path = str(harness.PARAMS_DIR / f"{ini}.ini")
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-I", "-c", WORKER, str(harness.tree_path(tree)), ini_path, mode, wantcls,
                            wantscal, wanttens, str(HERE.parent)], capture_output=True, text=True, timeout=900, cwd=OUT,
                           env=dict(os.environ, OMP_NUM_THREADS="1"), preexec_fn=lambda: os.nice(10))
    rc = p.returncode
    sig = signal.Signals(-rc).name if rc < 0 else ""
    (OUT / f"{name}.log").write_text(p.stdout + "\n--- stderr ---\n" + p.stderr)
    lines = [ln for ln in p.stdout.splitlines() if ln.strip()]
    err = [ln for ln in p.stderr.splitlines() if "Fortran runtime error" in ln or "At line" in ln or "Error" in ln]
    return {"rc": rc, "signal": sig, "stdout_tail": lines[-3:], "stderr_key": err[:4]}


def main() -> None:
    only = set(sys.argv[1:])
    res_file = harness.AUDIT_ROOT / "work" / "verify" / "getevol_indep.json"
    res = json.loads(res_file.read_text()) if res_file.exists() else {}
    for case in CASES:
        if only and case[0] not in only:
            continue
        res[case[0]] = run_case(*case)
        print(case[0], json.dumps(res[case[0]]), flush=True)
        res_file.write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

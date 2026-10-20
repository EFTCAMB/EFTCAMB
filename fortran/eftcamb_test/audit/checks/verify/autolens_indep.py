"""Independent check of finding auto-lens-accuracy-enables-nonlinear-for-eft.

1) Python API state (no CAMB computation) on trees A0 (audit), B0 (upstream CAMB v2), B1 (pre-rebase EFTCAMB):
   NonLinear / max_eta_k after read_ini + set_for_lmax(2500) and after camb.set_params(lmax=2500).
2) Ini path (Fortran reader, no set_for_lmax) through harness.run_model on A0: do_nonlinear=0 vs 2 at identical
   l_max_scalar / k_eta_max_scalar, for f(R) wCDM, f(R) LCDM, GR wCDM (EFTflag=0, w=-0.9), GR LCDM (EFTflag=0).
   Gives the size of the non-linear lensing correction with EFT on vs off.
3) Warning visibility: f(R) with do_nonlinear=2, EFTCAMB_feedback_level=0 vs 1, grep the log.
4) B1 (pre-rebase) opt-in non-linear lensing on f(R) wCDM: same NL/lin ratio as A0?

Run: /Users/marco/Software/miniforge3/bin/python autolens_indep.py
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REPO = HERE.parents[4]
PARAMS = REPO / "fortran" / "eftcamb_test" / "parameters"
SCR = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
TREES = SCR / "trees"
OUT = SCR / "work" / "verify" / "autolens"
OUT.mkdir(parents=True, exist_ok=True)

STATE_CODE = r"""
import sys, json, os
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, sys.argv[1])
import camb
from pathlib import Path
assert Path(camb.__file__).resolve().is_relative_to(Path(sys.argv[1]).resolve())
res = {"version": camb.__version__}
cp = camb.read_ini(sys.argv[2])
res["read_ini_NonLinear"] = str(cp.NonLinear)
cp.set_for_lmax(2500)
res["read_ini+set_for_lmax(2500)"] = [str(cp.NonLinear), float(cp.max_eta_k)]
try:
    cp2 = camb.set_params(H0=70, ombh2=0.0226, omch2=0.112, lmax=2500)
except Exception:
    cp2 = camb.CAMBparams(); cp2.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, YHe=0.24); cp2.set_for_lmax(2500)
    res["set_params_note"] = "camb.set_params failed (old code vs new numpy); used set_cosmology+set_for_lmax"
res["set_params(lmax=2500)"] = [str(cp2.NonLinear), float(cp2.max_eta_k)]
cp3 = camb.read_ini(sys.argv[2]); cp3.set_for_lmax(2500, lens_potential_accuracy=1)
res["set_for_lmax(2500,lpa=1)"] = [str(cp3.NonLinear), float(cp3.max_eta_k)]
print("JSON" + json.dumps(res))
"""


def api_state() -> dict:
    ini = OUT / "state.ini"
    ini.write_text(f"DEFAULT({PARAMS / '4_Des_F_R_wCDM_1.ini'})\nfeedback_level = 0\nEFTCAMB_write_background = F\n")
    out = {}
    for t in ("A0", "B0", "B1"):
        p = subprocess.run([sys.executable, "-c", STATE_CODE, str(TREES / t), str(ini)], capture_output=True,
                           text=True, cwd=OUT, timeout=300)
        line = [x for x in p.stdout.splitlines() if x.startswith("JSON")]
        out[t] = json.loads(line[0][4:]) if line else {"error": p.stderr[-800:]}
    return out


COMMON = {"l_max_scalar": "2700", "k_eta_max_scalar": "72000", "get_tensor_cls": "F", "get_transfer": "F",
          "feedback_level": "0", "EFTCAMB_feedback_level": "0"}
CASES = {
    "fR_wCDM": (PARAMS / "4_Des_F_R_wCDM_1.ini", {}),
    "fR_LCDM": (PARAMS / "4_Des_F_R_LCDM_1.ini", {}),
    "GR_wCDM": (PARAMS / "4_Des_F_R_wCDM_1.ini", {"EFTflag": "0", "w": "-0.9"}),
    "GR_LCDM": (PARAMS / "4_Des_F_R_wCDM_1.ini", {"EFTflag": "0"}),
}


def lensed(npz: Path) -> np.ndarray:
    return harness.load(npz)["cl_total"]


def ratio_stats(nl: np.ndarray, lin: np.ndarray, lmax: int = 2500) -> dict:
    res = {}
    for lab, i in (("TT", 0), ("EE", 1), ("BB", 2)):
        r = np.abs(nl[2 : lmax + 1, i] / lin[2 : lmax + 1, i] - 1)
        res[lab] = (float(r.max()), int(np.argmax(r) + 2), float(nl[lmax, i] / lin[lmax, i] - 1))
    return res


def ini_runs(tree: str = "A0", cases=None) -> dict:
    jobs = []
    for name, (ini, ov) in (cases or CASES).items():
        for nl in ("0", "2"):
            jobs.append((name, nl, ini, {**COMMON, **ov, "do_nonlinear": nl}))
    from concurrent.futures import ThreadPoolExecutor

    def one(j):
        name, nl, ini, ov = j
        npz = OUT / f"{tree}_{name}_nl{nl}.npz"
        return j, harness.run_model(tree, ini, npz, overrides=ov, threads=2)

    with ThreadPoolExecutor(4) as ex:
        sts = list(ex.map(one, jobs))
    res = {}
    for (name, nl, _, _), st in sts:
        if st["status"] != "ok":
            res[f"{name}_nl{nl}"] = st
    for name in (cases or CASES):
        try:
            res[name] = ratio_stats(lensed(OUT / f"{tree}_{name}_nl2.npz"), lensed(OUT / f"{tree}_{name}_nl0.npz"))
        except Exception as e:  # noqa: BLE001
            res[name] = repr(e)
    return res


def warning_check() -> dict:
    res = {}
    for fb in ("0", "1"):
        npz = OUT / f"A0_warn_fb{fb}.npz"
        ov = {**COMMON, "do_nonlinear": "2", "feedback_level": fb, "l_max_scalar": "1500",
              "k_eta_max_scalar": "4000"}
        st = harness.run_model("A0", PARAMS / "4_Des_F_R_wCDM_1.ini", npz, overrides=ov, threads=1)
        log = Path(st["log"]).read_text(errors="replace")
        res[f"EFTCAMB_feedback_level={fb}"] = {"status": st["status"],
                                               "warning_printed": "Nonlinearity has not been worked out" in log}
    return res


DICT_CODE = r"""
import sys, os
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, sys.argv[2]); import harness
sys.path.insert(0, sys.argv[1]); import camb
pars = camb.set_params(H0=70, ombh2=0.0226, omch2=0.112, lmax=2500, dark_energy_model="EFTCAMB")
eft = {"EFTflag": 3, "DesignerEFTmodel": 1, "EFTwDE": 1, "EFTw0": -0.9, "EFTB0": 1.0,
       "EFT_ghost_stability": True, "EFT_gradient_stability": True, "EFTCAMB_write_background": False}
pars.EFTCAMB.initialize_parameters(pars, eft, print_header=False)
print("NonLinear", pars.NonLinear, "EFTCAMB_feedback_level", pars.EFTCAMB.EFTCAMB_feedback_level, flush=True)
with harness.slot(1):
    camb.get_background(pars)
"""


def dict_warning() -> dict:
    p = subprocess.run([sys.executable, "-c", DICT_CODE, str(TREES / "A0"), str(HERE.parent)], capture_output=True,
                       text=True, cwd=OUT, timeout=600)
    txt = p.stdout + p.stderr
    return {"rc": p.returncode, "first_line": txt.splitlines()[0] if txt else "",
            "warning_printed": "Nonlinearity has not been worked out" in txt,
            "warning_line": [x for x in txt.splitlines() if "Warning: You have nonlinear" in x][:1]}


if __name__ == "__main__":
    what = sys.argv[1:] or ["state", "ini", "warn", "b1"]
    report = {}
    if "state" in what:
        report["api_state"] = api_state()
    if "ini" in what:
        report["ini_A0"] = ini_runs("A0")
    if "warn" in what:
        report["warning"] = warning_check()
    if "dict" in what:
        report["dict_warning"] = dict_warning()
    if "b1" in what:
        report["ini_B1"] = ini_runs("B1", {"fR_wCDM": CASES["fR_wCDM"], "GR_LCDM": CASES["GR_LCDM"]})
    print(json.dumps(report, indent=1))
    (OUT / f"indep_{'_'.join(what)}.json").write_text(json.dumps(report, indent=1))

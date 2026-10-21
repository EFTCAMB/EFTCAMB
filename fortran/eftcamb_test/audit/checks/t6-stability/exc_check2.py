"""(i) every instability reason -> clean Python exception; the same process then computes a stable model with results
identical to a fresh process.  Uses camb.get_results (full pipeline) in one subprocess per case."""
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness  # noqa: E402

TREE = str(harness.tree_path(sys.argv[1] if len(sys.argv) > 1 else "A0"))
P = harness.PARAMS_DIR
INIS = Path(__file__).resolve().parents[2] / "inis"
CASES = {
    "scalar ghost (PEFT gamma1=-1)": (P / "2_PEFT_gamma1_const_1.ini", {"EFTGamma10": "-1.0"}),
    "scalar gradient (RPH alpha_B=-0.3)": (INIS / "t6_RPH_KB_const.ini", {"RPHbraiding0": "-0.3"}),
    "tensor gradient (RPH alpha_T=-1.2)": (INIS / "t6_RPH_KB_const.ini", {"RPHtensor0": "-1.2", "RPHbraiding0": "0.0"}),
    "beyond GLPV with ghost check (PEFT gamma6)": (P / "2_PEFT_gamma6_const_1.ini", {"EFT_ghost_stability": "T"}),
    "mass (quintessence, mass check on)": (P / "5_quint_exponential.ini", {"EFT_mass_stability": "T"}),
    "additional prior (K-mouflage alphaU<0)": (P / "5_Kmouflage_1.ini", {"alphaU": "-1.0"}),
    "ghost_math sign change (window)": (INIS / "t6_PEFT_gamma1_taylor_window.ini",
                                        {"EFT_ghost_math_stability": "T", "EFT_ghost_stability": "F",
                                         "EFTGamma1a0": "0.5", "EFTGamma10": "-1e-3", "EFTGamma12": "1.0"}),
}
STABLE = P / "2_PEFT_gamma1_const_1.ini"
CODE = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
out = {}
cp = camb.read_ini(sys.argv[2], no_validate=True)
try:
    camb.get_results(cp)
    out["first"] = "no exception"
except Exception as e:
    out["first"] = type(e).__name__ + ": " + str(e).replace(chr(10), " | ")[:160]
cp2 = camb.read_ini(sys.argv[3], no_validate=True)
r2 = camb.get_results(cp2)
cl = r2.get_cmb_power_spectra(cp2, CMB_unit="muK")["total"][2:2000]
np.save(sys.argv[4], cl)
out["second_ok"] = True
print("T6EXC " + json.dumps(out))
"""
env = dict(os.environ, OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
OUTD = harness.AUDIT_ROOT / "work" / "t6-stability" / "exc_check"
OUTD.mkdir(parents=True, exist_ok=True)
SEL = [int(x) for x in sys.argv[2:]] if len(sys.argv) > 2 else list(range(len(CASES)))


def wrap(ini, ov, name):
    w = OUTD / f"{name}.ini"
    w.write_text("\n".join([f"DEFAULT({ini})", f"output_root = {OUTD}/{name}", "EFTCAMB_write_background = F",
                            "feedback_level = 0", *[f"{k} = {v}" for k, v in ov.items()]]) + "\n")
    return w


def run(args, timeout=420):
    try:
        with harness.slot(1):
            return subprocess.run([sys.executable, "-c", *args], env=env, cwd=OUTD, capture_output=True, text=True,
                                  timeout=timeout, preexec_fn=harness._nice)
    except subprocess.TimeoutExpired:
        return None


stable = wrap(STABLE, {}, "stable")
ref_file = OUTD / "ref.npy"
if not ref_file.exists():
    # fresh-process reference for the stable model
    ref_code = CODE.replace('cp = camb.read_ini(sys.argv[2], no_validate=True)', 'cp = None').replace(
        "    camb.get_results(cp)\n", "    pass\n")
    run([ref_code, TREE, str(stable), str(stable), str(ref_file)])
import numpy as np  # noqa: E402

ref = np.load(ref_file)
for i, (label, (ini, ov)) in enumerate(CASES.items()):
    if i not in SEL:
        continue
    w = wrap(ini, ov, f"case{i}")
    p = run([CODE, TREE, str(w), str(stable), str(OUTD / f"c{i}.npy")])
    if p is None:
        print(f"{label:45s} TIMEOUT (> 420 s): hang")
        continue
    res = [json.loads(l[6:]) for l in p.stdout.splitlines() if l.startswith("T6EXC ")]
    noise = [l for l in p.stdout.splitlines() if l.strip() and not l.startswith("T6EXC ")]
    if res:
        cl = np.load(OUTD / f"c{i}.npy")
        same = bool(np.array_equal(cl, ref))
        print(f"{label:45s} rc={p.returncode} first: {res[0]['first']}\n{'':45s} then stable model ok, Cl identical "
              f"to fresh process: {same}; stdout lines at feedback 0: {len(noise)} {noise[:2]}")
    else:
        print(f"{label:45s} rc={p.returncode} NO RESULT stderr: {p.stderr[-300:]!r}")

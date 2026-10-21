"""Independent check of finding minkowski-tachpos-lambdadotdot-zero (verifier).

Different path from the repro:
  * EFTLambdadotdot is rebuilt from the model's EFTLambda VALUES (not from EFTLambdadot), sampled exactly at the
    background interpolation-table nodes through get_eft_functions (11_EFTCAMB_python path, no stability code):
        Lambda = EFTLambda / a^2,  N = ln a,
        EFTLambdadot    = a^2 dLambda/dtau   = a^2 H Lambda_N                (convention check vs stored EFTLambdadot)
        EFTLambdadotdot = a^2 d2Lambda/dtau2 = a^2 (Hdot Lambda_N + H^2 Lambda_NN)
  * control models whose code sets EFTLambdadotdot analytically (5_scg full model, 2_PEFT_wDE_wCDM_1 designer) validate
    the reconstruction against the code's own analytic value.
  * the code's tachpos/ghostpos actually used by EFTTestStability are read through the camb ctypes structure
    (EFTCAMB_timestep_cache.from_buffer_copy) and the full Q1..Q9 / D1 formula of 06_abstract_EFTCAMB_model.f90 is
    re-implemented in numpy to attribute the code value to the EFTLambdadotdot it was fed.

Usage: [MINKTACH_MODELS=m1,m2] python minktach_indep.py [tree] [outjson]
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ROOT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
            "scratchpad/audit")
PAR = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")
TREE = sys.argv[1] if len(sys.argv) > 1 else str(ROOT / "trees/A0")
OUT = sys.argv[2] if len(sys.argv) > 2 else str(ROOT / "work/verify/minktach/indep_A0.json")
MODELS = ["5_quint_exponential", "5_quint_inverse_power", "5_quint_quadratic", "5_Kmouflage_1", "5_hdsk_TG_1",
          "5_scg", "2_PEFT_wDE_wCDM_1"]
if os.environ.get("MINKTACH_MODELS"):
    MODELS = os.environ["MINKTACH_MODELS"].split(",")

WORKER = r"""
import sys, json, ctypes
from ctypes import byref, c_double, c_int
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
from camb.baseconfig import camblib as lib
from camb.eftcamb import EFTCAMB_timestep_cache as C

ini = sys.argv[2]
cp = camb.read_ini(ini, no_validate=True)
cp.EFTCAMB.EFTCAMB_skip_stability = True
cp.EFTCAMB.EFT_mass_stability = False
res = camb.get_background(cp)
P = res.Params
eft = P.EFTCAMB
pars = eft.read_parameters()
designer = bool(eft.EFTCAMB_model_is_designer) if hasattr(eft, "EFTCAMB_model_is_designer") else None
# background table nodes (quintessence-type full models): x_i = x_ini + i dx
nnod = int(pars.get("model_background_num_points", 1100))
x0 = np.log(float(pars.get("model_background_a_ini", 1e-8)))
x1 = np.log(float(pars.get("model_background_a_final", 1.1)))
xn = x0 + np.arange(nnod) * (x1 - x0) / (nnod - 1)
sel = (xn > np.log(5e-4)) & (xn < 0.0)
a = np.exp(xn[sel])
fields, arr = eft.get_eft_functions(res, a)
g = lambda f: np.array(arr[f], dtype=float)
H, Hd = g("adotoa"), g("Hdot")
Lam = g("EFTLambda") / a**2
N = np.log(a)
LN = np.gradient(Lam, N, edge_order=2)
LNN = np.gradient(LN, N, edge_order=2)
Ld_val = a**2 * H * LN
Ldd_val = a**2 * (Hd * LN + H**2 * LNN)
# FD of stored EFTLambdadot (repro method) for comparison
Ld = g("EFTLambdadot")
Ldd_fdyp = H * np.gradient(Ld, N, edge_order=2) - 2 * H * Ld
Ldd_code_geteft = g("EFTLambdadotdot")

# code values used by EFTTestStability (minkowski on so the full branch is exercised; verdict per point)
eft.EFT_minkowski_limit = True
test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
buf = (ctypes.c_char * (1 << 18))()
need = ["adotoa", "Hdot", "Hdotdot", "Hdotdotdot", "EFTc", "EFTLambdadot", "EFTLambdadotdot", "EFTGamma2V",
        "EFTGamma2P", "EFTGamma2PP", "EFTGamma2PPP", "EFTGamma3V", "EFTGamma3P", "EFTGamma3PP", "EFTGamma3PPP",
        "EFTGamma3PPPP", "ghostpos", "tachpos"]
S = {f: np.empty_like(a) for f in need}
ok = np.empty(len(a), bool)
for i, x in enumerate(a):
    ok[i] = bool(test(byref(c_double(x)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf))
    c = C.from_buffer_copy(bytes(buf[:ctypes.sizeof(C)]))
    for f in need:
        S[f][i] = getattr(c, f)
h0 = float(P.EFTCAMB_parameter_cache.h0_Mpc)

def tach(s, Ldd):
    A = s["adotoa"]; Hd_ = s["Hdot"]; Hdd = s["Hdotdot"]; Hddd = s["Hdotdotdot"]; a2 = a**2
    D1 = (-2*s["EFTc"] + 3*a*A*h0*s["EFTGamma2V"] + (A**2 + 2*Hd_)*s["EFTGamma3V"] + 2*a*A**2*s["EFTGamma3P"])
    Q1 = (-A*s["EFTLambdadot"] + Ldd)*a2
    Q2 = h0*(Hdd - A*Hd_)*a2*s["EFTGamma2P"]/2
    Q3 = 3*h0*A*Hd_*a**3*s["EFTGamma2PP"]/2
    Q4 = h0*A**3*a2**2*s["EFTGamma2PPP"]/2
    Q5 = (-Hddd - 12*A**4 + 18*A**2*Hd_ - 2*Hd_**2 + A*Hdd)*s["EFTGamma3V"]
    Q6 = (a/2)*(-Hddd + 12*A**4 - 5*Hd_**2 - 5*A*Hdd)*s["EFTGamma3P"]
    Q7 = -(a2/2)*(3*Hd_**2 + 8*A**2*Hd_ + 4*Hdd*A)*s["EFTGamma3PP"]
    Q8 = -a**3*A**2*(3*Hd_ + A**2)*s["EFTGamma3PPP"]
    Q9 = -a2**2*A**4*s["EFTGamma3PPPP"]/2
    return D1, -(Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9)/D1

D1_py, tach_py_code = tach(S, S["EFTLambdadot"]*0 + S["EFTLambdadotdot"])
_, tach_py_val = tach(S, Ldd_val)
_, tach_py_fdyp = tach(S, Ldd_fdyp)
s = slice(3, -3)
late = s
def rel(x, y):
    return float(np.max(np.abs(x[s] - y[s])) / np.max(np.abs(y[s])))
fail = lambda D1, t: (-D1 <= 0) | (-t < 0)
i5 = int(np.argmin(abs(a - 0.5)))
out = dict(
    model_name=None, n=int(len(a[s])), a_range=[float(a[s][0]), float(a[s][-1])], dx=float(xn[1]-xn[0]),
    conv_check_rel_Ld_val_vs_stored=rel(Ld_val, Ld),
    max_abs_Ldd_code_geteft=float(np.abs(Ldd_code_geteft[s]).max()),
    max_abs_Ldd_code_stab=float(np.abs(S["EFTLambdadotdot"][s]).max()),
    rel_Ldd_val_vs_code_geteft=(rel(Ldd_val, Ldd_code_geteft) if np.abs(Ldd_code_geteft[s]).max() > 0 else None),
    rel_Ldd_val_vs_fdyp=rel(Ldd_val, Ldd_fdyp),
    rel_ghostpos_py_vs_code=rel(D1_py, S["ghostpos"]),
    rel_tachpos_py_vs_code=rel(tach_py_code, S["tachpos"]),
    max_g2g3=float(np.abs(S["EFTGamma2V"]).max() + np.abs(S["EFTGamma3V"]).max()),
    max_abs_g3PP_g2PP=float(np.abs(S["EFTGamma3PP"]).max() + np.abs(S["EFTGamma2PP"]).max()),
    per_point_unstable_code=int((~ok)[s].sum()),
    fail_code=int(fail(S["ghostpos"], S["tachpos"])[s].sum()),
    fail_with_Ldd_val=int(fail(S["ghostpos"], tach_py_val)[s].sum()),
    fail_with_Ldd_fdyp=int(fail(S["ghostpos"], tach_py_fdyp)[s].sum()),
    a05=dict(a=float(a[i5]), Ldd_stab=float(S["EFTLambdadotdot"][i5]), Ldd_geteft=float(Ldd_code_geteft[i5]),
             Ldd_val=float(Ldd_val[i5]), Ldd_fdyp=float(Ldd_fdyp[i5]), tachpos_code=float(S["tachpos"][i5]),
             tachpos_Ldd_val=float(tach_py_val[i5]), ghostpos=float(S["ghostpos"][i5])),
)
print("RESULT " + json.dumps(out), flush=True)
"""


def main() -> None:
    rows = {}
    with tempfile.TemporaryDirectory() as d:
        for m in MODELS:
            ini = Path(d) / f"{m}.ini"
            ini.write_text(f"DEFAULT({PAR}/{m}.ini)\nfeedback_level = 0\nEFTCAMB_write_background = F\n"
                           f"output_root = {d}/out\nEFTCAMB_stability_time = 1.e-3\nEFTCAMB_feedback_level = 0\n")
            env = dict(os.environ, OMP_NUM_THREADS="1")
            env.pop("PYTHONPATH", None)
            with harness.slot(1):
                p = subprocess.run([sys.executable, "-c", WORKER, TREE, str(ini)], capture_output=True, text=True,
                                   cwd=d, env=env, preexec_fn=lambda: os.nice(10), timeout=1200)
            res = [json.loads(line[7:]) for line in p.stdout.splitlines() if line.startswith("RESULT ")]
            if not res:
                rows[m] = {"error": (p.stdout[-1500:] + p.stderr[-1500:])}
                print(m, "ERROR", rows[m]["error"][-800:])
                continue
            rows[m] = r = res[0]
            print(f"== {m}: n={r['n']} a in {r['a_range']} dx={r['dx']:.4f}")
            print(f"   convention: |a^2 H Lambda_N - EFTLambdadot|/max = {r['conv_check_rel_Ld_val_vs_stored']:.2e}")
            print(f"   max|Ldd| code(get_eft_functions) = {r['max_abs_Ldd_code_geteft']:.3e}; "
                  f"code(at positivity, stability cache) = {r['max_abs_Ldd_code_stab']:.3e}")
            print(f"   Ldd(values) vs code analytic: {r['rel_Ldd_val_vs_code_geteft']}; vs FD(EFTLambdadot): "
                  f"{r['rel_Ldd_val_vs_fdyp']:.2e}")
            print(f"   python formula vs code: ghostpos {r['rel_ghostpos_py_vs_code']:.1e}, tachpos "
                  f"{r['rel_tachpos_py_vs_code']:.1e}; max|g2|+|g3| = {r['max_g2g3']:.2e}; "
                  f"max|g2PP|+|g3PP| = {r['max_abs_g3PP_g2PP']:.2e}")
            print(f"   Minkowski fails: code {r['fail_code']}/{r['n']} (per-point unstable {r['per_point_unstable_code']}),"
                  f" with Ldd(values) {r['fail_with_Ldd_val']}/{r['n']}, with FD(EFTLambdadot) "
                  f"{r['fail_with_Ldd_fdyp']}/{r['n']}")
            print(f"   a~0.5: {r['a05']}")
    Path(OUT).parent.mkdir(parents=True, exist_ok=True)
    Path(OUT).write_text(json.dumps({"tree": TREE, "rows": rows}, indent=1))
    print("tree:", TREE, "->", OUT)


if __name__ == "__main__":
    main()

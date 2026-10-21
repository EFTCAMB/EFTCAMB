"""t3-gr: checks with tree R8 = A0 + pid_d2 (EFTpiDfunction = D1 + k^2 D2) + getevol_scalars (get_evolution forces
WantScalars while setting up the k modes).
 1) impact of the pi-dot IC fix on spectra (A0 vs R8, same settings) for models with D1/D2 /= 0;
 2) get_scale_evolution after a transfer-only run on R8 (subprocess; A0 segfaults) and agreement with the evolution
    obtained after a run with C_l.
"""
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

from common import HERE, PARAMS, RUNS, WORK, fmt, harness

SUB = RUNS / "c6"
R8 = WORK / "R8"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "1500", "num_redshiftwindows": "0"}
MODELS = ["3_RPH_MK_const_1", "5_Horava_eta_1", "2_PEFT_Omega_const_1", "4_Des_F_R_LCDM_1", "5_hdsk_TG_1"]

WORKER = r'''
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, sys.argv[1])
import camb, numpy as np
cp = camb.read_ini(sys.argv[2], no_validate=True)
cp.WantCls = sys.argv[3] == "cls"
cp.WantScalars = cp.WantCls
cp.WantTensors = False
cp.InitPower.r = 0.0
cp.DoLensing = False
cp.max_l = 300
cp.WantTransfer = True
cp.SourceWindows = []  # transfer-only + source windows is an upstream crash (t10a)
r = camb.get_results(cp)
a = np.logspace(-3, 0, 40)
names, arr = cp.EFTCAMB.get_scale_evolution(r, [0.05, 0.17], a)
out = {f: np.asarray(arr[f]).tolist() for f in ("pi", "pidot", "clxc", "z") if f in arr.dtype.names}
print("JSON" + json.dumps(out))
'''


def run_cmp(item):
    tree, m = item
    out = SUB / f"{'R8' if tree == R8 else 'A0'}_{m}.npz"
    if not out.with_suffix(".json").exists():
        harness.run_model(tree, PARAMS / f"{m}.ini", out, overrides=FAST, threads=1, timeout=540)
    return tree, m, out


def evol(tree, mode):
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-I", "-c", WORKER, str(tree), str(PARAMS / "3_RPH_MK_const_1.ini"), mode],
                           capture_output=True, text=True, timeout=600, cwd=SUB)
    js = [ln[4:] for ln in p.stdout.splitlines() if ln.startswith("JSON")]
    return p.returncode, (json.loads(js[0]) if js else None)


def main() -> None:
    SUB.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(run_cmp, [(t, m) for m in MODELS for t in ("A0", R8)]))
    outs = {(("R8" if t == R8 else "A0"), m): o for t, m, o in res}
    print("1) pi-dot IC fix (R8) vs A0, same settings:")
    for m in MODELS:
        a, b = outs[("A0", m)], outs[("R8", m)]
        sa = json.loads(a.with_suffix(".json").read_text())["status"]
        sb = json.loads(b.with_suffix(".json").read_text())["status"]
        print(f"   {m:24s} {sa}/{sb}", fmt(harness.compare(b, a, lmax=1450)) if sa == sb == "ok" else "")
    print("2) get_scale_evolution after transfer-only run (3_RPH_MK_const_1, k=0.05,0.17, a=1e-3..1):")
    for tree in ("A0", R8):
        tp = harness.tree_path(tree)
        rc_t, ev_t = evol(tp, "transfer")
        rc_c, ev_c = evol(tp, "cls")
        line = f"   {('R8' if tree == R8 else 'A0'):3s} transfer-only rc={rc_t} ({'ok' if ev_t else 'crashed'}); with C_l rc={rc_c}"
        if ev_t and ev_c:
            import numpy as np
            d = max(float(np.max(np.abs(np.array(ev_t[f]) - np.array(ev_c[f])) / (np.max(np.abs(ev_c[f])) + 1e-300)))
                    for f in ev_c)
            line += f"; max rel diff transfer-only vs C_l run over {sorted(ev_c)} = {d:.1e}"
        print(line)


if __name__ == "__main__":
    main()

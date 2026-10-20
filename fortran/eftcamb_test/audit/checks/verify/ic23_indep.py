"""Independent check of ic23 finding: harness runs, IC1 vs IC2/3 at late (1e-2) and early (1e-8) turn-on, A0 and B1."""
import sys, json
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
import numpy as np
P = harness.REPO / "fortran/eftcamb_test/parameters"
OUT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/ic23"
jobs = []
for tree in ("A0", "B1"):
    for ton in ("1e-2", "1e-8"):
        for ic in (1, 2):
            jobs.append((tree, "2_PEFT_Omega_const_2", ton, ic))
for ton in ("1e-2", "1e-8"):
    for ic in (1, 3):
        jobs.append(("A0", "2_PEFT_gamma3_const_2", ton, ic))
def go(j):
    tree, m, ton, ic = j
    out = f"{OUT}/{tree}_{m}_on{ton}_ic{ic}.npz"
    st = harness.run_model(tree, P / f"{m}.ini", out,
                           overrides={"EFT_IC_type": str(ic), "EFTCAMB_turn_on_time": ton}, evol_k=[0.01, 0.1, 1.0])
    return j, out, st["status"], st.get("log_tail", "")[-200:]
with ThreadPoolExecutor(4) as ex:
    res = list(ex.map(go, jobs))
st = {r[0]: (r[1], r[2], r[3]) for r in res}
for j, (o, s, t) in st.items():
    if s != "ok": print(j, s, t)
def cmp(a, b, label):
    if st[a][1] != "ok" or st[b][1] != "ok":
        print(label, "skip", st[a][1], st[b][1]); return
    A, B = harness.load(st[a][0]), harness.load(st[b][0])
    c = harness.compare(A, B)
    ea = np.array(A["evol_a"]); i = int(np.argmax(ea > float(a[2]) * 1.03))
    pa, pb = A["evol_pi"][:, i], B["evol_pi"][:, i]
    print(f"{label}: worst_cl={c['worst_cl']:.3e} worst_pk={c['worst_pk']:.3e} s8 {float(A['sigma8'][-1]) if np.ndim(A['sigma8']) else float(A['sigma8']):.5f}->"
          f"{float(B['sigma8'][-1]) if np.ndim(B['sigma8']) else float(B['sigma8']):.5f}  pi@a={ea[i]:.3g}: {pa} vs {pb}")
for tree in ("A0", "B1"):
    for ton in ("1e-2", "1e-8"):
        cmp((tree, "2_PEFT_Omega_const_2", ton, 1), (tree, "2_PEFT_Omega_const_2", ton, 2), f"{tree} Omega_const_2 on={ton} IC1 vs IC2")
for ton in ("1e-2", "1e-8"):
    cmp(("A0", "2_PEFT_gamma3_const_2", ton, 1), ("A0", "2_PEFT_gamma3_const_2", ton, 3), f"A0 gamma3_const_2 on={ton} IC1 vs IC3")

import sys, os
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness, numpy as np
OUT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/ic23"
def f(tree, m, ton, ic): return f"{OUT}/{tree}_{m}_on{ton}_ic{ic}.npz"
def cmp(a, b, label, ton):
    if not (os.path.exists(a) and os.path.exists(b)): print(label, "missing", os.path.exists(a), os.path.exists(b)); return
    A, B = harness.load(a), harness.load(b); c = harness.compare(A, B)
    ea = A["evol_a"][0]; i = int(np.argmax(ea > max(float(ton), 1e-6) * 1.03))
    print(f"{label}: worst_cl={c['worst_cl']:.3e} worst_pk={c['worst_pk']:.3e} s8 {float(A['sigma8'][0]):.5f}->{float(B['sigma8'][0]):.5f}"
          f" pi(k=.01,.1,1)@a={ea[i]:.3g}: {np.array2string(A['evol_pi'][:, i], precision=3)} vs {np.array2string(B['evol_pi'][:, i], precision=3)}")
for tree in ("A0", "B1"):
    for ton in ("1e-2", "1e-8"):
        cmp(f(tree, "2_PEFT_Omega_const_2", ton, 1), f(tree, "2_PEFT_Omega_const_2", ton, 2), f"{tree} Omega_const_2 on={ton} IC1 vs IC2", ton)
cmp(f("A0", "2_PEFT_gamma3_const_2", "1e-8", 1), f("A0", "2_PEFT_gamma3_const_2", "1e-8", 3), "A0 gamma3_const_2 on=1e-8 IC1 vs IC3", "1e-8")
cmp(f("A0", "2_PEFT_Omega_const_2", "1e-2", 2), f("B1", "2_PEFT_Omega_const_2", "1e-2", 2), "IC2 on=1e-2 A0 vs B1", "1e-2")

import os, sys, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
from scan_all import sampler_grid
INSTR = str(t6lib.WORK / "instr")
log = t6lib.WORK / "instr_rph_mb.log"
os.environ["T6_STAB_LOG"] = str(log)
ini = harness.PARAMS_DIR / "3_RPH_MB_OmegaDE_1.ini"
lin, log1, _ = sampler_grid()
r = t6lib.probe(INSTR, ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0"}, a_grid=log1,
                feedback=0, cleanup_each=False, check=True)
print("per-point (ctypes) unstable on log1:", (~r["verdict"]).sum(), " full_check:", r["full_check"])
rows = [l.split() for l in log.read_text().splitlines()]
loops = np.array([int(x[0]) for x in rows]); A = np.array([float(x[1]) for x in rows]); S = np.array([x[2] == "T" for x in rows])
K = np.array([float(x[3]) for x in rows])
for L in (1, 2, 3):
    m = loops == L
    print(f"loop {L}: n={m.sum()} n_fail={(~S[m]).sum()} a_range=[{A[m].min() if m.any() else None}, {A[m].max() if m.any() else None}] unique_a={len(np.unique(A[m]))}")
m2 = loops == 2
if m2.any():
    # compare loop-2 points with python's log1 grid
    n = min(m2.sum(), len(log1))
    print("max rel diff code-vs-python loop-2 a:", np.max(np.abs(A[m2][:n] - log1[:n]) / log1[:n]))
    kin_ct = r["EFT_kinetic"][:n]
    print("loop-2 kinetic: code<0:", (K[m2][:n] < 0).sum(), " ctypes<0:", (kin_ct < 0).sum(),
          " identical values:", int(np.sum(K[m2][:n] == kin_ct)))
    i = np.where(kin_ct < 0)[0][:3]
    print(" sample a, kin(code loop), kin(ctypes):", [(A[m2][j], K[m2][j], kin_ct[j]) for j in i])

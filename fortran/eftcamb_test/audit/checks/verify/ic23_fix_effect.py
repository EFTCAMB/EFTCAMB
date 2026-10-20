"""Does the candidate cache fill (EFT_T5_ICFIX=1 in t5-robust/A0_icdbg) change IC3 outputs / remove NaNs, and is the
ICFIX=0 path of that tree identical to frozen A0?"""
import sys, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import harness
S = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
tree = S + "/work/t5-robust/A0_icdbg"
out = pathlib.Path(S + "/work/verify/ic23"); out.mkdir(parents=True, exist_ok=True)
P = harness.REPO / "fortran/eftcamb_test/parameters"
for m in sys.argv[1:]:
    for fix in ("0", "1"):
        for thr in (1, 2):
            f = out / f"icdbg_{m}_on1e-8_ic3_fix{fix}_t{thr}.npz"
            st = harness.run_model(tree, str(P / f"{m}.ini"), str(f),
                                   overrides={"EFT_IC_type": 3, "EFTCAMB_turn_on_time": "1e-8"}, threads=thr,
                                   extra_env={"EFT_T5_ICFIX": fix, "EFT_T5_DEBUG": "0"})
            if st["status"] != "ok":
                print(m, fix, thr, st["status"], st.get("log_tail", "")[-3:]); continue
            d = harness.load(str(f)); c = d["cl_total"]
            print(m, "fix", fix, "threads", thr, "nNaN", int(np.isnan(c).sum()), "TT[10,100,1000]", c[[10, 100, 1000], 0],
                  "sigma8", d["sigma8"], flush=True)

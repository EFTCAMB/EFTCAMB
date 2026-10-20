"""Verifier: run 3_RPH_MB_OmegaDE_1 through harness with stability_time overrides, feedback 3, report verdict and messages."""
import sys
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
W = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/rphkin"
ini = harness.PARAMS_DIR / "3_RPH_MB_OmegaDE_1.ini"
for tree in ("A0", "B1"):
    for st in ("0.1", "1e-5", "6e-6", "1e-6", "1e-10"):
        r = harness.run_model(tree, ini, f"{W}/{tree}_st{st}.npz", overrides={"EFTCAMB_stability_time": st,
                              "EFTCAMB_feedback_level": "3"}, threads=1)
        log = open(r["log"]).read()
        lines = [l for l in log.splitlines() if "instab" in l.lower() or "unstable" in l.lower() or "ghost" in l.lower()]
        print(tree, st, r["status"], lines[:4])

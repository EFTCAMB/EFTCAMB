import sys, numpy as np
sys.path.insert(0, ".")
import t6lib
FIX = str(t6lib.WORK / "fixy")
ini = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/inis/t6_PEFT_gamma1_taylor_window.ini"
for tree in ("A0", FIX):
    r = t6lib.probe(tree, ini, overrides={"feedback_level": "0"})
    print(tree[-6:], "get_background with stability on:", r["status"], r.get("bg_error"))

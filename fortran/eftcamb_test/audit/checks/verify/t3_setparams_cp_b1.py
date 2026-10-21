"""Check whether pre-rebase EFTCAMB (B1) has the same set_params(cp=...) reset (EFTflag only, no CAMB run).

Run under a slot: python run_in_slot.py 1 python t3_setparams_cp_b1.py <tree>
"""

import sys

sys.path.insert(0, sys.argv[1])
import camb  # noqa: E402

eft = {"EFTflag": 3, "DesignerEFTmodel": 1, "EFTwDE": 0, "EFTB0": 0.1}
p = camb.set_params(H0=67.5, ombh2=0.0224, omch2=0.12, YHe=0.245, **eft)
print(sys.argv[1].rstrip("/").split("/")[-1], "before:", p.EFTCAMB.EFTflag, end="  ")
p = camb.set_params(cp=p, WantTensors=False)
print("after set_params(cp=p, WantTensors=False):", p.EFTCAMB.EFTflag)

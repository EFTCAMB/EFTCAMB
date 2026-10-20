"""Child for an lldb backtrace: build a designer f(R) EFT model, deepcopy the params, then call `action`.

Usage: python pickle_segfault_lldb_child.py TREE background|num_params [INI]
"""

import copy
import sys

sys.path.insert(0, sys.argv[1])
import camb  # noqa: E402

if len(sys.argv) > 3:  # read_ini route (Fortran ini reader: usable on -fcheck trees, avoids F-001)
    p = camb.read_ini(sys.argv[3], no_validate=True)
else:
    p = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965)
    p.EFTCAMB.initialize_parameters(p, {"EFTflag": 3, "DesignerEFTmodel": 1, "EFTwDE": 0, "EFTB0": 1.0},
                                    print_header=False)
q = copy.deepcopy(p)
print("after deepcopy EFTflag =", q.EFTCAMB.EFTflag, flush=True)
if sys.argv[2] == "background":
    camb.get_background(q)
else:
    q.EFTCAMB.num_params()
print("no crash", flush=True)

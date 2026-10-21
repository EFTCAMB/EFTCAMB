"""Child for the lldb probe: get_background for a sequence of models in one process (release tree).

UNINIT_ALL_BG=1 uses get_background for every model (same call depth); otherwise
get_results/get_background alternate. A model spec is NAME or NAME@key=value@key=value (ini overrides).
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, sys.argv[1])
import camb  # noqa: E402

PARAMS = Path(__file__).resolve().parents[3] / "parameters"
work = Path(sys.argv[2])
for i, spec in enumerate(sys.argv[3:]):
    model, *ov = spec.split("@")
    ini = work / f"{model}_{i}.ini"
    ini.write_text("feedback_level = 0\nEFTCAMB_feedback_level = 0\nEFTCAMB_write_background = F\n"
                   "get_transfer = F\nget_tensor_cls = F\ndo_lensing = F\n" + "".join(f"{o}\n" for o in ov) +
                   f"DEFAULT({PARAMS / (model + '.ini')})\n")
    print(f"MODEL {spec}", flush=True)
    cp = camb.read_ini(str(ini))
    try:
        camb.get_results(cp) if (i % 2 and os.environ.get("UNINIT_ALL_BG") != "1") else camb.get_background(cp)
        print("  accepted", flush=True)
    except Exception as e:  # noqa: BLE001
        print(f"  ERROR {e}", flush=True)

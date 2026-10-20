"""t3-gr (b): approach to GR of Pure EFT constant Gamma1 (2_PEFT_gamma1_const_1, EFT_mass_stability=T).

The Gamma1 -> 0 scan hits 'Mass instability' at Gamma1 = 1e-6 although Gamma1 = 1, 1e-2 are stable and 1e-9 is
switched to GR by the RGR logic. Print the reported mu_2/H^2 (or 'stable') versus Gamma1 to get the scaling.
"""

import re
import sys
from concurrent.futures import ThreadPoolExecutor

from common import PARAMS, harness, RUNS

INI = PARAMS / "2_PEFT_gamma1_const_1.ini"
VALS = [1e-1, 1e-2, 3e-3, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 3e-8]
OV = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "50", "get_transfer": "F", "num_redshiftwindows": "0"}


def one(v):
    out = RUNS / "b4" / f"G1_{v:g}.npz"
    st = harness.run_model("A0", INI, out, overrides={**OV, "EFTGamma10": repr(v)}, threads=1, timeout=900)
    log = out.with_suffix(".log").read_text(errors="replace")
    mu = re.findall(r"Mass instability\. (mu_\d/H\^2) = *([-+0-9.E]+)", log)
    at = re.findall(r"Instability detected at a = *([-+0-9.E]+)", log)
    rgr = re.findall(r"Return to GR time: *([-+0-9.E]+)", log)
    return v, st["status"], mu[:1], at[:1], rgr[:1]


with ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 3) as ex:
    for v, s, mu, at, rgr in ex.map(one, VALS):
        print(f"Gamma1={v:8.1e} status={s:12s} RGR={rgr} {mu} at a={at}")

"""Extra checks for t3-gamma1-mass-instability-near-gr.

a) pi(k, a) for constant Gamma1 = 1e-3 with the mass test off: a genuine mu^2 ~ -500 H^2 would make pi grow like
   exp(22 N); measure max |pi| and |pidot| at k = 0.01, 0.1 1/Mpc.
b) scope probe: same float cancellation for another near-GR model? designer f(R) on LCDM with small B0 and
   EFT_mass_stability = T.
"""

import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "gamma1_mu2"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "50", "get_transfer": "F",
        "num_redshiftwindows": "0"}


def job(spec):
    tag, ini, ov, evk = spec
    out = OUT / f"{tag}.npz"
    st = harness.run_model("A0", harness.PARAMS_DIR / ini, out, overrides={**FAST, **ov}, evol_k=evk, threads=1,
                           timeout=900)
    log = out.with_suffix(".log").read_text(errors="replace")
    mu = re.findall(r"Mass instability\. mu_(\d)/H\^2 = *([-+0-9.E]+)", log)
    at = re.findall(r"Instability detected at a = *([-+0-9.E]+)", log)
    return tag, st["status"], mu[:1], at[:1], out


def main() -> None:
    specs = [("evol_G0.001_nomass", "2_PEFT_gamma1_const_1.ini",
              {"EFTGamma10": "0.001", "EFT_mass_stability": "F"}, [0.01, 0.1])]
    for b0 in ("1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7"):
        specs.append((f"fR_B0_{b0}_mass", "4_Des_F_R_LCDM_1.ini", {"EFTB0": b0, "EFT_mass_stability": "T"}, None))
    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(job, specs))
    for tag, s, mu, at, out in res:
        print(f"{tag:22s} {s:12s} mass-msg={mu} at={at}")
        if tag.startswith("evol") and s == "ok":
            d = np.load(out, allow_pickle=True)
            for f in ("pi", "pidot"):
                key = "evol_" + f
                if key in d:
                    arr = d[key]
                    print(f"   max|{f}| over k={list(d['evol_k'])}, a in [1e-6,1]: {np.nanmax(np.abs(arr)):.3e}")
            for f in ("clxc",):
                if "evol_" + f in d:
                    print(f"   max|{f}|: {np.nanmax(np.abs(d['evol_' + f])):.3e}")


if __name__ == "__main__":
    main()

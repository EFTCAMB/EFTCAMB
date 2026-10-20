"""Repro (t3-gr): Pure EFT constant Gamma1 is rejected as mass-unstable as Gamma1 -> 0 (approach to GR blocked).

Model 2_PEFT_gamma1_const_1.ini (EFTflag=1, all EFT functions zero except a constant Gamma1, w=-1, LCDM background,
EFT_mass_stability=T). The pi action then only has 2 M_2^4 (pi_dot + delta N)^2: it is shift symmetric, so no
mass term should appear, and as Gamma1 -> 0 the model must tend to LCDM. EFTCAMB's mass eigenvalue mu_2
(06_abstract_EFTCAMB_model.f90, EFT_mu2, divided by F1^2 F2^4 F3^2) instead grows roughly like 1/Gamma1^2, always at
the first sampled time after the turn-on a = 0.01; Gamma1 = 1e-2 is accepted, 1e-3 ... 3e-8 are rejected, and
Gamma1 <= 1e-8 is accepted again because the RGR logic switches the model to GR.
Expected: accepted (or a Gamma1-independent mu_2/H^2) for all Gamma1. Observed: printed.
Run: /Users/marco/Software/miniforge3/bin/python t3-gamma1-mass-instability-near-gr.py
"""

import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

INI = harness.PARAMS_DIR / "2_PEFT_gamma1_const_1.ini"
OUT = harness.AUDIT_ROOT / "work" / "t3-gr" / "repro_gamma1_mass"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "50", "get_transfer": "F",
        "num_redshiftwindows": "0"}
VALS = [1e-2, 1e-3, 1e-5, 1e-7, 1e-9]


def job(v):
    out = OUT / f"G1_{v:g}.npz"
    st = harness.run_model("A0", INI, out, overrides={**FAST, "EFTGamma10": repr(v)}, threads=1, timeout=900)
    log = out.with_suffix(".log").read_text(errors="replace")
    mu = re.findall(r"Mass instability\. mu_\d/H\^2 = *([-+0-9.E]+)", log)
    at = re.findall(r"Instability detected at a = *([-+0-9.E]+)", log)
    rgr = re.findall(r"Return to GR time: *([-+0-9.E]+)", log)
    return v, st["status"], mu[:1], at[:1], rgr[:1]


def main() -> None:
    print(f"{'Gamma1':>8s} {'status':12s} {'RGR time':>11s} {'mu_2/H^2':>12s} {'at a':>11s}")
    with ThreadPoolExecutor(3) as ex:
        for v, s, mu, at, rgr in ex.map(job, VALS):
            print(f"{v:8.0e} {s:12s} {rgr[0] if rgr else '-':>11s} {mu[0] if mu else '-':>12s} {at[0] if at else '-':>11s}")
    print("expected: all ok (constant Gamma1 on LCDM has no mass term); observed: rejected for 1e-3 ... 1e-7")


if __name__ == "__main__":
    main()

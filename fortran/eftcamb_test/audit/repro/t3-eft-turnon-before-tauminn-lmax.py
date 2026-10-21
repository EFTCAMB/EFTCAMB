"""Repro (t3-gr): an EFT model whose perturbations turn on early (5_ADE_1: EFTCAMB_turn_on_time = 1e-10) is refused
with a hard error as soon as l_max_scalar is lowered to 2000, because the thermal-history start depends on kmax.

equations.f90 (initial, "EFTCAMB MOD START: compute the time at which turn on EFT", ~l.2186) raises
GlobalError('EFTCAMB starts before thermo tauminn ...') when tau(EFTCAMB_turn_on_time) < ThermoData%tauminn.
cmbmain sets taumin = GetTauStart(maxq) (maxq from the CMB/transfer k range) and ThermoData%tauminn = 0.95 taumin
without taking the EFT turn-on into account, so a smaller l_max (smaller k_max, later tau start) makes a valid
model fail: tau(a=1e-10) = 4.6e-5 Mpc vs tauminn = 5.0e-5 at l_max_scalar = 2000 (ok at 2500, 3500).
Expected: the run succeeds for any l_max (thermal history started at min(GetTauStart(maxq), tau_turn_on)).
Run: /Users/marco/Software/miniforge3/bin/python t3-eft-turnon-before-tauminn-lmax.py   (2 CAMB runs, ~1 min)
"""

import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

INI = harness.PARAMS_DIR / "5_ADE_1.ini"
OUT = harness.AUDIT_ROOT / "work" / "t3-gr" / "repro_tauminn"


def job(lmax):
    out = OUT / f"ADE_lmax{lmax}.npz"
    st = harness.run_model("A0", INI, out, overrides={"l_max_scalar": str(lmax), "get_tensor_cls": "F"}, threads=1,
                           timeout=900)
    msg = re.findall(r"(EFTCAMB starts before thermo tauminn[^\n]*)", out.with_suffix(".log").read_text(errors="replace"))
    return lmax, st["status"], msg[:1]


def main() -> None:
    with ThreadPoolExecutor(2) as ex:
        for lmax, status, msg in ex.map(job, [2500, 2000]):
            print(f"l_max_scalar={lmax}: {status} {msg}")
    print("expected: ok for both; observed: l_max_scalar=2000 refused")


if __name__ == "__main__":
    main()
